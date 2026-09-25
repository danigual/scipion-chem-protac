# **************************************************************************
# *
# * Name:     test of protocol_protac_model.py
# *
# * Authors: Daniel Gutiérrez (Internship student, github user: danigual)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Tests for ProtPROTACModel and its driver.

The driver runs under Python 2 in production; here its helpers are exercised under the
test runner's Python 3, with PROTAC-Model's own modules replaced by stubs. The only test
that needs the real tools (TestPROTACModelExample1) skips itself when they are missing,
since the test runner cannot skip on a *_HOME variable.
"""

import argparse
import contextlib
import importlib.util
import io
import os
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock

from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from protac import Plugin
from protac.protocols import ProtPROTACModel

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.join(os.path.dirname(_TESTS_DIR), 'scripts', 'run_protac_model.py')
_TOY_PDB = os.path.join(_TESTS_DIR, 'data', 'transplant', 'toy_source.pdb')


def _write(path, content=''):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    return path


def _read(path):
    with open(path) as f:
        return f.read()


@contextlib.contextmanager
def _chdir(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


class _FakePipe:
    """ What os.popen() returns: read() and a close() status in os.wait() format. """
    def __init__(self, output, status=None):
        self._output, self._status = output, status

    def read(self):
        return self._output

    def close(self):
        return self._status


class _FakeOs:
    def __init__(self, output='', status=None):
        self.commands = []
        self._output, self._status = output, status
        self.sep = os.sep

    def popen(self, cmd, *args, **kwargs):
        self.commands.append(cmd)
        return _FakePipe(self._output, self._status)


class TestPROTACModelDriver(unittest.TestCase):
    """ The driver's own code: the Vina fixes, the fail-safe pose filtering and the
    staging each phase does around PROTAC-Model's functions. """

    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('run_protac_model_undertest', _SCRIPT)
        cls.driver = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.driver)

    def setUp(self):
        self.workDir = tempfile.mkdtemp(prefix='protac_model_driver_')
        self.addCleanup(shutil.rmtree, self.workDir, True)
        # Every test starts with PROTAC-Model not loaded.
        self.driver.pre = self.driver.fro = self.driver.ros = None
        self.driver._warningCounts.clear()

    def _pdbqt(self, *coords, extra=''):
        lines = [extra]
        for x, y, z in coords:
            lines.append('HETATM    1  C   UNL     1    %8.3f%8.3f%8.3f  1.00  0.00' % (x, y, z))
        return _write(os.path.join(self.workDir, 'lig.pdbqt'), '\n'.join(lines) + '\n')

    # ------------------------------ Vina box ------------------------------
    def test_readPdbqtCoords_skips_other_and_malformed_records(self):
        path = self._pdbqt((1, 2, 3), extra='REMARK  VINA RESULT\n'
                                             'ATOM      2  C   UNL     1      xxxxx')
        self.assertEqual(self.driver._readPdbqtCoords(path), [(1.0, 2.0, 3.0)])

    def test_vinaBoxFlags_pads_the_bounding_box(self):
        path = self._pdbqt((0, 0, 0), (10, 2, 0))
        # Centre of the bounding box; half-extent + 4 A padding, rounded up, doubled.
        self.assertEqual(self.driver._vinaBoxFlags(path),
                         '--center_x 5.000 --center_y 1.000 --center_z 0.000 '
                         '--size_x 18.000 --size_y 10.000 --size_z 8.000')

    def test_vinaBoxFlags_without_atoms(self):
        self.assertIsNone(self.driver._vinaBoxFlags(self._pdbqt()))
        self.assertIsNone(self.driver._vinaBoxFlags(os.path.join(self.workDir, 'no.pdbqt')))

    def test_isVinaScoreOnly(self):
        scoreOnly = 'vina --score_only --receptor r.pdbqt --ligand l.pdbqt | grep Affinity'
        self.assertTrue(self.driver._isVinaScoreOnly(scoreOnly))
        for cmd in ('vina --receptor r --ligand l', scoreOnly + ' --center_x 1',
                    scoreOnly + ' --autobox', scoreOnly + ' --maps m'):
            self.assertFalse(self.driver._isVinaScoreOnly(cmd), cmd)

    def test_vina_score_regex_reads_both_output_formats(self):
        regex = self.driver._VINA_SCORE_RE
        self.assertEqual(regex.search('Affinity: -7.5 (kcal/mol)\n').group(1), '-7.5')
        self.assertEqual(regex.search(
            'Intramolecular energy: -1.2\n'
            'Estimated Free Energy of Binding   : -8.123 (kcal/mol) [=(1)+(2)+(3)+(4)]\n'
            '(1) Final Intermolecular Energy    : -9.0 (kcal/mol)\n').group(1), '-8.123')
        self.assertIsNone(regex.search('    Ligand - Receptor : -9.0\n'))

    # ------------------------------ Vina rerun ------------------------------
    def _scoreOnly(self, fakeOs, ligand=None):
        ligand = ligand or self._pdbqt((0, 0, 0), (10, 2, 0))
        cmd = f'vina --score_only --receptor r.pdbqt --ligand {ligand} | grep Affinity | cut -f2'
        return self.driver._runVinaScoreOnly(fakeOs, cmd, (), {})

    def test_runVinaScoreOnly_adds_a_box_and_parses_the_score(self):
        fakeOs = _FakeOs('Estimated Free Energy of Binding   : -8.1 (kcal/mol)\n')
        self.assertEqual(self._scoreOnly(fakeOs), '-8.1\n')
        # The grep|cut tail is dropped, and stderr is kept for the warning.
        self.assertEqual(len(fakeOs.commands), 1)
        self.assertNotIn('|', fakeOs.commands[0])
        self.assertIn('--center_x 5.000', fakeOs.commands[0])
        self.assertTrue(fakeOs.commands[0].endswith('2>&1'))

    def test_runVinaScoreOnly_scores_zero_on_failure(self):
        with mock.patch('sys.stdout', io.StringIO()):
            self.assertEqual(self._scoreOnly(_FakeOs('ERROR: outside the grid box')), '0\n')
            self.assertEqual(self._scoreOnly(_FakeOs('x'), ligand=self._pdbqt()), '0\n')

    def test_runVinaScoreOnly_stops_when_vina_cannot_run(self):
        with self.assertRaises(self.driver._FatalError):
            self._scoreOnly(_FakeOs('vina: not found', status=127 << 8))

    def test_runVinaScoreOnly_ignores_an_unexpected_command(self):
        self.assertIsNone(self.driver._runVinaScoreOnly(
            _FakeOs(), 'vina --score_only --ligand lig.pdb', (), {}))

    def test_guardedOs_replaces_empty_output_with_zero(self):
        guarded = self.driver._GuardedOs(_FakeOs(''))
        with mock.patch('sys.stdout', io.StringIO()):
            self.assertEqual(guarded.popen('grep -c x f').read(), '0\n')
        self.assertEqual(self.driver._GuardedOs(_FakeOs('3\n')).popen('wc').read(), '3\n')
        # Everything else is the real module.
        self.assertEqual(guarded.sep, os.sep)

    # ------------------------------ pose filtering ------------------------------
    def test_failSafeFiltering_skips_a_failing_pose(self):
        def filtering(paraList):
            raise IndexError('no conformer')

        module = types.SimpleNamespace(filtering=filtering)
        wrapped = self.driver._makeFailSafeFiltering(module, 0)
        # pool.map() pickles it by name.
        self.assertEqual(wrapped.__name__, 'filtering')
        with _chdir(self.workDir), mock.patch('sys.stdout', io.StringIO()) as out:
            self.assertIsNone(wrapped([17, 'x']))
            self.driver._reportPoseFailures()
            self.assertEqual(_read(self.driver.POSE_FAILURES_LOG), '17\n')
        self.assertIn('1 pose(s) skipped', out.getvalue())

    def test_failSafeFiltering_lets_fatal_errors_through(self):
        def filtering(paraList):
            raise self.driver._FatalError('vina missing')

        wrapped = self.driver._makeFailSafeFiltering(types.SimpleNamespace(filtering=filtering), 0)
        with self.assertRaises(self.driver._FatalError):
            wrapped([1])

    def test_installRuntimeFixes_without_rosetta(self):
        self.driver.pre = types.SimpleNamespace(os=os)
        self.driver.fro = types.SimpleNamespace(filtering=lambda p: 'ok')
        with mock.patch.dict(os.environ, {'LC_ALL': 'es_ES.UTF-8'}):
            self.driver.installRuntimeFixes()
            self.driver.installRuntimeFixes()
            # awk under a decimal-comma locale would drop every conformer.
            self.assertEqual(os.environ['LC_ALL'], 'C')
        self.assertIsInstance(self.driver.pre.os, self.driver._GuardedOs)
        # Not wrapped twice.
        self.assertIs(self.driver.pre.os._os, os)
        self.assertEqual(self.driver.fro.filtering([1]), 'ok')

    def test_loadPipeline_imports_rosetta_only_to_refine(self):
        home = os.path.join(self.workDir, 'protac-model')
        _write(os.path.join(home, 'utils', '__init__.py'))
        _write(os.path.join(home, 'utils', 'preprocess.py'), 'NAME = "pre"\n')
        _write(os.path.join(home, 'utils', 'frodock.py'), 'NAME = "fro"\n')
        # Like the real one, it needs ROSETTA as soon as it is imported.
        _write(os.path.join(home, 'utils', 'rosetta.py'),
               'import os\nROSETTA = os.environ["ROSETTA"]\n')
        pathBefore = list(sys.path)
        modules = {k: v for k, v in sys.modules.items() if k == 'utils' or k.startswith('utils.')}
        try:
            with mock.patch.dict(os.environ, {'PROTAC_MODEL_HOME': home}), \
                    mock.patch.dict(sys.modules):
                for name in modules:
                    del sys.modules[name]
                os.environ.pop('ROSETTA', None)
                self.driver.loadPipeline(withRosetta=False)
                self.assertEqual((self.driver.pre.NAME, self.driver.fro.NAME), ('pre', 'fro'))
                self.assertIsNone(self.driver.ros)
                with self.assertRaises(KeyError):
                    self.driver.loadPipeline(withRosetta=True)
        finally:
            sys.path[:] = pathBefore

    # ------------------------------ phases ------------------------------
    def test_runFrodock_stages_its_inputs(self):
        calls = []
        self.driver.pre = types.SimpleNamespace(
            alter_pro_chain=lambda *a: calls.append(('alter',) + a))
        self.driver.fro = types.SimpleNamespace(frodock=lambda site: calls.append(('frodock', site)))
        frodockHome = os.path.join(self.workDir, 'frodock')
        _write(os.path.join(frodockHome, 'bin', 'soap.bin'), 'soap')
        e3 = [_write(os.path.join(self.workDir, f'e3_{i}.sdf'), f'e3 {i}') for i in (1, 2)]
        run = os.path.join(self.workDir, 'run')
        _write(os.path.join(run, 'frodock_score.txt'), 'stale')
        args = argparse.Namespace(receptor='/in/rec.pdb', target='/in/tgt.pdb',
                                  smiles=' CCO ', site='1.0,2.0,3.0', e3lig1=e3[0], e3lig2=e3[1])
        with mock.patch.dict(os.environ, {'FRODOCK': frodockHome}), _chdir(run):
            self.driver.runFrodock(args)
            self.assertEqual(_read('protac.smi'), 'CCO\n')
            self.assertEqual((_read('rec_lig_1.sdf'), _read('rec_lig_2.sdf')), ('e3 1', 'e3 2'))
            self.assertEqual(_read('soap.bin'), 'soap')
            # frodock() appends to it.
            self.assertFalse(os.path.exists('frodock_score.txt'))
        self.assertEqual(calls, [('alter', '/in/rec.pdb', '/in/tgt.pdb', 'receptor.pdb',
                                  'target.pdb'), ('frodock', '1.0,2.0,3.0')])

    def test_runFilter_clears_appended_files_first(self):
        seen = []

        def filterFrodock(cpu, ligLocateNum, targetSmi, recSmi):
            seen.append(((cpu, ligLocateNum, targetSmi, recSmi),
                         sorted(os.path.relpath(os.path.join(d, f))
                                for d, _, files in os.walk('.') for f in files)))

        self.driver.fro = types.SimpleNamespace(filter_frodock=filterFrodock)
        for name in ('results_voromqa', 'addH_log', 'vina/score_all_top1', 'vina/score_filter',
                     'rec_lig_1/vina/score_filter', 'keep.pdb'):
            _write(os.path.join(self.workDir, name), 'old')
        args = argparse.Namespace(cpu=4, ligLocateNum=1, targetSmi='none', recSmi='CCO')
        with _chdir(self.workDir):
            self.driver.runFilter(args)
        self.assertEqual(seen, [((4, 1, 'none', 'CCO'),
                                 ['keep.pdb', os.path.join('rec_lig_1', 'vina', 'score_filter')])])

    def test_runRefine_copies_what_rosetta_needs(self):
        seen = []
        self.driver.ros = types.SimpleNamespace(rosetta=lambda *a: seen.append(sorted(os.listdir('.'))))
        for name in ('rec_lig.sdf', 'target_lig.sdf', 'protac.smi', 'rec_lig_1.sdf', 'rec_lig_2.sdf'):
            _write(os.path.join(self.workDir, 'frodock', name), name)
        rosettaDir = os.path.join(self.workDir, 'rosetta')
        os.makedirs(rosettaDir)
        with _chdir(rosettaDir):
            self.driver.runRefine(argparse.Namespace(cpu=12, ligLocateNum=2, targetSmi='none',
                                                     recSmi='none'))
        self.assertEqual(seen, [['protac.smi', 'rec_lig.sdf', 'rec_lig_1.sdf', 'rec_lig_2.sdf',
                                 'target_lig.sdf']])


class TestPROTACModelProtocol(BaseTest):
    """ The protocol's own code, with steps called directly on a saved protocol. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        protImport = cls.newProtocol(ProtImportPdb, inputPdbData=1, pdbFile=_TOY_PDB)
        cls.launchProtocol(protImport)
        cls.structure = protImport.outputPdb

        dataDir = protImport._getExtraPath()
        e3Set = SetOfSmallMolecules().create(dataDir, suffix='e3')
        for name in ('e3_1', 'e3_2', 'dup', 'dup'):
            e3Set.append(SmallMolecule(smallMolFilename=_write(
                os.path.join(dataDir, f'{name}.sdf'), f'{name}\n'), molName=name))
        protImport._defineOutputs(e3Ligands=e3Set)
        protImport._store()
        cls.e3Ligands = protImport.e3Ligands

    def _newProt(self, e3=None, **overrides):
        args = {'inputReceptor': self.structure, 'inputTarget': self.structure,
                'siteCoords': '1.5,-2,3', 'protacSmiles': 'CCO'}
        if e3 is not None:
            args.update(e3Ligands=self.e3Ligands, e3Ligand1Name=e3[0], e3Ligand2Name=e3[1])
        args.update(overrides)
        return self.newProtocol(ProtPROTACModel, **args)

    def _newSavedProt(self, **overrides):
        prot = self._newProt(**overrides)
        self.saveProtocol(prot)
        prot.makeWorkingDir()
        return prot

    # ------------------------------ helpers ------------------------------
    def test_getSiteCoords(self):
        cases = {'1.5,-2,3': (1.5, -2.0, 3.0), ' 1 , 2.5 ,3e1 ': (1.0, 2.5, 30.0),
                 '1,2': None, '1,2,3,4': None, 'a,b,c': None, '': None, None: None}
        for value, expected in cases.items():
            self.assertEqual(self._newProt(siteCoords=value)._getSiteCoords(), expected, value)

    def test_getSiteArg_has_no_spaces(self):
        self.assertEqual(self._newProt(siteCoords=' 1 , 2.5 ,-3 ')._getSiteArg(),
                         '1.0000,2.5000,-3.0000')
        self.assertEqual(self._newProt(siteCoords='1,2')._getSiteArg(), '1,2')

    def test_e3_conformers_only_count_as_a_pair(self):
        prot = self._newProt()
        self.assertEqual((prot._getLigLocateNum(), prot._getE3LigandArgs()), (1, [0, '', '']))
        self.assertEqual(self._newProt(e3=('e3_1', ''))._getLigLocateNum(), 1)
        prot = self._newProt(e3=(' e3_1 ', 'e3_2'))
        self.assertEqual(prot._getLigLocateNum(), 2)
        self.assertEqual(prot._getE3LigandArgs(),
                         [self.e3Ligands.getObjId(), 'e3_1', 'e3_2'])

    def test_convertInputStep_writes_the_e3_conformers(self):
        prot = self._newSavedProt(e3=('e3_2', 'e3_1'))
        with mock.patch.object(Plugin, 'prepareFrodockBinDir') as prepareFrodock:
            prot.convertInputStep(_TOY_PDB, 'LIG', _TOY_PDB, '', *prot._getE3LigandArgs())
        prepareFrodock.assert_called_once_with(prot._getFrodockBinDir())
        self.assertEqual([_read(f) for f in prot._getE3LigandFiles()], ['e3_2\n', 'e3_1\n'])
        self.assertTrue(os.path.exists(prot._getReceptorFile()))
        self.assertTrue(os.path.exists(prot._getTargetFile()))

    def test_getSmiArg(self):
        for value, expected in ((None, 'none'), ('', 'none'), (' CCO ', 'CCO')):
            prot = self._newProt(targetLigandSmiles=value)
            self.assertEqual(ProtPROTACModel._getSmiArg(prot.targetLigandSmiles), expected)

    def test_getSmiArg_blank_means_none(self):
        prot = self._newProt(targetLigandSmiles='   ')
        self.assertEqual(ProtPROTACModel._getSmiArg(prot.targetLigandSmiles), 'none')

    def test_cleanReceptorOrTarget_keeps_only_the_named_warhead(self):
        pdb = ('ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
               'HETATM    2  C1  LIG A 101       1.000   0.000   0.000  1.00  0.00           C\n'
               'HETATM    3  S   SO4 A 102       2.000   0.000   0.000  1.00  0.00           S\n'
               'HETATM    4  O   HOH A 103       3.000   0.000   0.000  1.00  0.00           O\n'
               'END\n')
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        inFile = _write(os.path.join(tmp, 'in.pdb'), pdb)

        def resnames(ligandName):
            outFile = os.path.join(tmp, 'out.pdb')
            ProtPROTACModel._cleanReceptorOrTarget(inFile, outFile, ligandName)
            return {l[17:20] for l in _read(outFile).splitlines()
                    if l.startswith(('ATOM', 'HETATM'))}

        # Lower case on purpose: residue names are upper-cased.
        self.assertEqual(resnames(' lig '), {'ALA', 'LIG'})
        # Without a name only water goes.
        self.assertEqual(resnames(''), {'ALA', 'LIG', 'SO4'})

    # ------------------------------ validation ------------------------------
    def test_validate_rejects_a_malformed_site(self):
        for value in ('1,2', 'a,b,c'):
            errors = self._newProt(siteCoords=value)._validate()
            self.assertTrue(any('3 comma-separated' in e for e in errors), errors)
        self.assertFalse(any('3 comma-separated' in e for e in self._newProt()._validate()))

    def test_validate_handles_an_empty_site(self):
        # An empty form field arrives as ''.
        errors = self._newProt(siteCoords='').validate()
        self.assertTrue(any('Receptor interface site' in e for e in errors), errors)

    def test_validate_refinement_needs_more_than_10_threads(self):
        errors = self._newProt(doRefine=True, numberOfThreads=10)._validate()
        self.assertTrue(any('more than 10 threads' in e for e in errors))
        errors = self._newProt(doRefine=True, numberOfThreads=11)._validate()
        self.assertFalse(any('more than 10 threads' in e for e in errors))

    def test_validate_e3_conformer_names(self):
        cases = {('e3_1', ''): 'needs both conformer names',
                 ('e3_1', 'e3_1'): 'must be different',
                 ('e3_1', 'nope'): 'no molecule named "nope"',
                 ('dup', 'e3_1'): '2 molecules are named "dup"'}
        for names, fragment in cases.items():
            errors = self._newProt(e3=names)._validate()
            self.assertTrue(any(fragment in e for e in errors), (names, errors))
        errors = self._newProt(e3=('e3_1', 'e3_2'))._validate()
        self.assertFalse(any('E3' in e for e in errors), errors)

    def test_validate_reports_each_error_once(self):
        errors = self._newProt().validate()
        self.assertEqual(len(errors), len(set(errors)), errors)

    # ------------------------------ step arguments ------------------------------
    def _stepArgs(self, **overrides):
        prot = self._newProt(**overrides)
        prot._insertAllSteps()
        return [(step.funcName.get(), step.argsStr.get()) for step in prot._steps]

    def test_changing_an_input_reruns_its_step_and_the_later_ones(self):
        base = self._stepArgs()
        self.assertEqual([name for name, _ in base], ['convertInputStep', 'frodockStep',
                                                      'filterPosesStep', 'createOutputStep'])
        for field, value, firstChanged in (('receptorLigandName', 'LIG', 0),
                                           ('siteCoords', '1,2,3', 1),
                                           ('protacSmiles', 'CCN', 1),
                                           ('targetLigandSmiles', 'C', 2)):
            changed = self._stepArgs(**{field: value})
            for i, (old, new) in enumerate(zip(base, changed)):
                if i < firstChanged:
                    self.assertEqual(old, new, (field, old[0]))
                else:
                    self.assertNotEqual(old, new, (field, old[0]))

    def test_threads_only_rerun_the_refinement(self):
        # Without refinement they only size the process pool, so nothing is redone.
        self.assertEqual(self._stepArgs(numberOfThreads=4), self._stepArgs(numberOfThreads=8))
        base = self._stepArgs(doRefine=True, numberOfThreads=12)
        changed = self._stepArgs(doRefine=True, numberOfThreads=16)
        self.assertEqual(base[:3], changed[:3])
        for old, new in zip(base[3:], changed[3:]):
            self.assertNotEqual(old, new, old[0])

    def test_adding_e3_conformers_reruns_everything(self):
        base, changed = self._stepArgs(), self._stepArgs(e3=('e3_1', 'e3_2'))
        for old, new in zip(base, changed):
            self.assertNotEqual(old, new, old[0])

    def test_refinement_adds_its_step(self):
        names = [name for name, _ in self._stepArgs(doRefine=True, numberOfThreads=12)]
        self.assertEqual(names, ['convertInputStep', 'frodockStep', 'filterPosesStep',
                                 'refineStep', 'createOutputStep'])

    # ------------------------------ output ------------------------------
    def _results(self, prot, lines, models, refine=False):
        folder = 'rosetta_results' if refine else 'frodock_results'
        resultsDir = prot._getExtraPath(folder, 'all')
        name = 'results_rosetta.txt' if refine else 'results_frodock.txt'
        if lines is not None:
            _write(os.path.join(resultsDir, name), lines)
        for poseId in models:
            _write(os.path.join(resultsDir, f'model_merge_{poseId}.pdb'), _read(_TOY_PDB))

    def _outputScores(self, prot):
        return {os.path.basename(m.getFileName()): m._score.get()
                for m in prot.outputTernaryModels}

    def test_createOutputStep_one_item_per_listed_pose(self):
        prot = self._newSavedProt()
        # A pose without its model file is skipped.
        self._results(prot, '3 -120.5\n\n7 -98.25\n9 -50\n', models=(3, 7))
        prot.createOutputStep('[]')
        self.assertEqual(self._outputScores(prot), {'model_merge_3.pdb': -120.5,
                                                    'model_merge_7.pdb': -98.25})

    def test_createOutputStep_reads_the_refined_results(self):
        prot = self._newSavedProt(doRefine=True, numberOfThreads=12)
        self._results(prot, '1 -1.0\n', models=(1,))
        self._results(prot, '2 -2.0\n', models=(2,), refine=True)
        prot.createOutputStep('[]')
        self.assertEqual(self._outputScores(prot), {'model_merge_2.pdb': -2.0})

    def test_createOutputStep_fails_without_usable_poses(self):
        for lines, models in ((None, ()), ('', ()), ('1 -5.0\n', ())):
            prot = self._newSavedProt()
            self._results(prot, lines, models)
            with self.assertRaises(RuntimeError):
                prot.createOutputStep('[]')

    def test_createOutputStep_skips_a_malformed_line(self):
        prot = self._newSavedProt()
        self._results(prot, '1 -5.0\n2\n', models=(1, 2))
        prot.createOutputStep('[]')
        self.assertEqual(self._outputScores(prot), {'model_merge_1.pdb': -5.0})


class TestPROTACModelExample1(BaseTest):
    """ example_1 of the PROTAC-Model repository (VHL + SMARCA2 with its PROTAC), run
    end to end. It ships with the PROTAC-Model checkout, so it is read from there. Takes
    hours; skips itself when any of the tools is not installed. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def _exampleDir(self):
        try:
            Plugin.getProtacModelEnviron()
            Plugin.getProtacModelPython()
        except FileNotFoundError as e:
            self.skipTest(f'PROTAC-Model or one of its tools is not installed: {e}')
        exampleDir = os.path.join(Plugin.getProtacModelScript(), 'example_1')
        if not os.path.isdir(exampleDir):
            self.skipTest(f'{exampleDir} not found')
        return exampleDir

    def test_example_1(self):
        exampleDir = self._exampleDir()
        structures = {}
        for name in ('receptor', 'target'):
            protImport = self.newProtocol(ProtImportPdb, inputPdbData=1,
                                          pdbFile=os.path.join(exampleDir, f'{name}.pdb'))
            self.launchProtocol(protImport)
            structures[name] = protImport.outputPdb

        prot = self.newProtocol(
            ProtPROTACModel, inputReceptor=structures['receptor'], receptorLigandName='FWZ',
            inputTarget=structures['target'], targetLigandName='FWZ',
            siteCoords=_read(os.path.join(exampleDir, 'site_info.txt')).strip(),
            protacSmiles=_read(os.path.join(exampleDir, 'protac.smi')).strip(),
            targetLigandSmiles=_read(os.path.join(exampleDir, 'target_lig.smi')).strip(),
            numberOfThreads=8)
        self.launchProtocol(prot)

        models = getattr(prot, 'outputTernaryModels', None)
        self.assertIsNotNone(models, 'No outputTernaryModels was generated')
        self.assertGreater(len(models), 0)
        for model in models:
            self.assertTrue(os.path.exists(model.getFileName()))
            self.assertIsNotNone(model._score.get())
