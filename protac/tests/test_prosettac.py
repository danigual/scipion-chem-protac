# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Name:     test of protocol_prosettac.py
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
Tests for ProtPRosettaC and its driver.

Everything here runs without PatchDock/Rosetta/PRosettaC installed. PRosettaC's own
modules (rosetta, protac_lib, clustering) are replaced by stubs, and PatchDock by three
small shell scripts, so what is tested is the driver's glue around them: file layout,
command lines, argument translation and failure handling.

What is NOT covered: an end-to-end run. It needs a real pair of structures with their
bound warheads plus a matching PROTAC SMILES, and hours of PatchDock + Rosetta.
"""

import argparse
import contextlib
import glob
import importlib.util
import json
import os
import shutil
import stat
import sys
import tempfile
import types
import unittest
from unittest import mock

from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from protac.protocols import ProtPRosettaC

_TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.join(os.path.dirname(_TESTS_DIR), 'scripts', 'run_prosettac.py')
_PATCHDOCK_OUTPUT = os.path.join(_TESTS_DIR, 'data', 'prosettac', 'Patchdock_output')
_TOY_PDB = os.path.join(_TESTS_DIR, 'data', 'transplant', 'toy_source.pdb')


def _loadDriver(scriptsFol):
    """ Imports run_prosettac.py as a module. It is shipped as plugin data, not as part of
    an importable package (same as pwchem's own scripts/), so it goes through importlib.
    SCRIPTS_FOL must exist first: the driver reads it at import time and prepends it to
    sys.path. Pointing it at an empty directory keeps PRosettaC's own rosetta.py/utils.py
    from shadowing the same-named modules of the scipion-chem-rosetta plugin. """
    os.environ['SCRIPTS_FOL'] = scriptsFol
    spec = importlib.util.spec_from_file_location('run_prosettac_undertest', _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write(path, content=''):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    return path


def _read(path):
    with open(path) as f:
        return f.read()


def _writeExecutable(path, content):
    _write(path, content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@contextlib.contextmanager
def _chdir(path):
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


class TestPRosettaCDriver(unittest.TestCase):
    """ The driver's own code - the part that is ours rather than PRosettaC's. """

    @classmethod
    def setUpClass(cls):
        cls._tmpRoot = tempfile.mkdtemp(prefix='prosettac_driver_')
        cls._scriptsFol = os.path.join(cls._tmpRoot, 'scripts_fol')
        os.makedirs(cls._scriptsFol)
        cls._pathBefore = list(sys.path)
        cls._scriptsFolBefore = os.environ.get('SCRIPTS_FOL')
        cls.driver = _loadDriver(cls._scriptsFol)

    @classmethod
    def tearDownClass(cls):
        # Both point into the temporary tree removed below; leaving either behind would
        # follow this process into every later test.
        sys.path[:] = cls._pathBefore
        if cls._scriptsFolBefore is None:
            os.environ.pop('SCRIPTS_FOL', None)
        else:
            os.environ['SCRIPTS_FOL'] = cls._scriptsFolBefore
        shutil.rmtree(cls._tmpRoot, ignore_errors=True)

    def setUp(self):
        self.workDir = tempfile.mkdtemp(prefix='prosettac_case_', dir=self._tmpRoot)

    def _writeState(self, state):
        _write(os.path.join(self.workDir, 'state.json'), json.dumps(state))

    def _readState(self):
        return json.loads(_read(os.path.join(self.workDir, 'state.json')))

    # ------------------------------ small helpers ------------------------------
    def test_dockingSuffix_matches_the_originals_own_format(self):
        # clustering.py parses the two numbers back out, so the leading zeros of the
        # docking index must be dropped.
        self.assertEqual(self.driver._dockingSuffix('pd.7_docking_0003.pdb'), '7_3')
        self.assertEqual(self.driver._dockingSuffix('pd.12_docking_0050.pdb'), '12_50')

    def test_stageSimpleName_copies_under_a_dot_free_name(self):
        src = _write(os.path.join(self.workDir, 'a.b.c', 'input.struct.pdb'), 'ATOM\n')
        with _chdir(self.workDir):
            self.assertEqual(self.driver._stageSimpleName(src, 'Struct0.pdb'), 'Struct0.pdb')
            self.assertEqual(_read('Struct0.pdb'), 'ATOM\n')

    def test_catFiles_concatenates_in_order(self):
        parts = [_write(os.path.join(self.workDir, f'part{i}.pdb'), content)
                 for i, content in enumerate(('first\n', 'second\n'))]
        dest = os.path.join(self.workDir, 'joined.pdb')
        self.driver._catFiles(parts, dest)
        # Order matters: it fixes the atom numbering PatchDock later sees.
        self.assertEqual(_read(dest), 'first\nsecond\n')

    def test_removeStale_only_removes_matching_files(self):
        for name in ('local.fasc', 'pd.1_docking_0001.pdb', 'keepme.pdb'):
            _write(os.path.join(self.workDir, name))
        self.driver._removeStale(self.workDir, 'local.fasc', '*_docking_????.pdb')
        self.assertEqual(os.listdir(self.workDir), ['keepme.pdb'])
        # A missing directory is not an error: phases clear leftovers before anything exists.
        self.driver._removeStale(os.path.join(self.workDir, 'nope'), '*.pdb')

    def test_require_rejects_missing_and_empty_files(self):
        empty = _write(os.path.join(self.workDir, 'empty.pdb'))
        for path in (empty, os.path.join(self.workDir, 'missing.pdb')):
            with self.assertRaises(RuntimeError):
                self.driver._require(path, 'test step')
        self.driver._require(_write(os.path.join(self.workDir, 'ok.pdb'), 'ATOM\n'), 'x')

    # ------------------------------ score.sc headers ------------------------------
    def test_dropRepeatedScoreHeaders_keeps_one_header_block_first(self):
        # clustering.py skips exactly the first two lines, so they must be the headers.
        scoreFile = _write(os.path.join(self.workDir, 'score.sc'),
                           'SEQUENCE: \n'
                           'SCORE: total_score fa_atr description\n'
                           'SCORE: -10.0 -1.0 combined_1_1_0001\n'
                           'SEQUENCE: \n'
                           'SCORE: total_score fa_atr description\n'
                           'SCORE: -20.0 -2.0 combined_2_1_0001\n')
        self.driver._dropRepeatedScoreHeaders(scoreFile)
        self.assertEqual(_read(scoreFile).splitlines(),
                         ['SEQUENCE: ',
                          'SCORE: total_score fa_atr description',
                          'SCORE: -10.0 -1.0 combined_1_1_0001',
                          'SCORE: -20.0 -2.0 combined_2_1_0001'])

    def test_dropRepeatedScoreHeaders_ignores_a_missing_file(self):
        self.driver._dropRepeatedScoreHeaders(os.path.join(self.workDir, 'score.sc'))

    # ------------------------------ PatchDock params ------------------------------
    def test_patchdockParams_puts_the_restraint_inline(self):
        text = ('receptorPdb Init0.pdb\n'
                '#distanceConstraints rec_atom_index lig_atom_index dist_thr\n'
                '#distanceConstraintsFile file_name\n'
                'clusterParams 0.1 4 2.0 4.0\n')
        out = self.driver._patchdockParams(text, [16, 25], 14, 2.0).splitlines()
        self.assertIn('receptorPdb Init0.pdb', out)
        self.assertIn('#distanceConstraintsFile file_name', out)
        self.assertIn('clusterParams 0.1 4 2.0 2.0', out)
        self.assertEqual(out[-1], 'distanceConstraints 16 25 14')

    def test_patchdockParams_disables_active_restraints(self):
        out = self.driver._patchdockParams('distanceConstraintsFile cst\n'
                                           'distanceConstraints 1 2 3\n'
                                           'clusterParams 0.1 4 2.0 4.0\n',
                                           [1, 2], 10.5, 2.0).splitlines()
        self.assertIn('#distanceConstraintsFile cst', out)
        self.assertIn('#distanceConstraints 1 2 3', out)
        # Only ours stays active.
        self.assertEqual([l for l in out if l.startswith('distanceConstraints')],
                         ['distanceConstraints 1 2 10.5'])

    def test_patchdockParams_rejects_an_unknown_clusterParams_line(self):
        with self.assertRaises(RuntimeError):
            self.driver._patchdockParams('clusterParams 0.1 4 2.0\n', [1, 2], 10, 2.0)

    # ------------------------------ PatchDock output ------------------------------
    def test_selectPatchdockSolutions_on_real_output(self):
        # Real rows, distances 7.41 / 7.93 / 13.54. The parameter header PatchDock writes
        # above them is replaced by a couple of lines that are not solution rows.
        output = _write(os.path.join(self.workDir, 'Patchdock_output'),
                        'receptorPdb Init0.pdb\n'
                        ' # | score | dist. || Ligand Transformation\n'
                        + _read(_PATCHDOCK_OUTPUT))
        select = self.driver._selectPatchdockSolutions
        self.assertEqual(select(output, 7.5, 500), [2, 3])
        self.assertEqual(select(output, 7.0, 2), [1, 2])
        self.assertEqual(select(output, 20, 5), [])

    def test_selectPatchdockSolutions_rejects_an_unreadable_distance(self):
        output = _write(os.path.join(self.workDir, 'Patchdock_output'),
                        '   1 | 13534 | -3.80 | n/a || 0 0 0 0 0 0\n')
        with self.assertRaises(RuntimeError):
            self.driver._selectPatchdockSolutions(output, 0, 5)

    def test_contiguousRanges_groups_consecutive_numbers(self):
        self.assertEqual(self.driver._contiguousRanges([1, 2, 3, 7, 8, 10]),
                         [(1, 3), (7, 8), (10, 10)])
        self.assertEqual(self.driver._contiguousRanges([]), [])

    def _fakePatchdock(self):
        """ buildParams.pl / patch_dock.Linux / transOutput.pl stand-ins. patch_dock
        returns the real output fixture; transOutput writes one model per solution, with
        the second warhead on chain X as PatchDock leaves it. """
        home = os.path.join(self.workDir, 'patchdock')
        _writeExecutable(os.path.join(home, 'buildParams.pl'),
                         '#!/bin/bash\n'
                         'printf "receptorPdb %s\\nligandPdb %s\\n'
                         '#distanceConstraintsFile file_name\\n'
                         'clusterParams 0.1 4 2.0 4.0\\n" "$1" "$2" > params.txt\n')
        _writeExecutable(os.path.join(home, 'patch_dock.Linux'),
                         '#!/bin/bash\n'
                         f'cp "$1" used_params.txt\ncp "{_PATCHDOCK_OUTPUT}" "$2"\n')
        _writeExecutable(os.path.join(home, 'transOutput.pl'),
                         '#!/bin/bash\n'
                         'echo "$2 $3" >> transOutput_calls.txt\n'
                         'for i in $(seq "$2" "$3"); do\n'
                         '  printf "REMARK solution %s\\nHETATM    1  C1  PT1 X   1\\n" "$i"'
                         ' > "$1.$i.pdb"\n'
                         'done\n')
        return home

    def _runPatchdock(self, minValue, maxValue=14.0, globalResults=500):
        self._writeState({'initFiles': ['Init0.pdb', 'Init1.pdb'], 'anchors': [15, 24],
                          'minValue': minValue, 'maxValue': maxValue})
        args = argparse.Namespace(globalResults=globalResults, threshold=2.0)
        with mock.patch.dict(os.environ, {'PATCHDOCK': self._fakePatchdock()}), \
                _chdir(self.workDir):
            self.driver.runPatchdock(args)

    def test_runPatchdock_keeps_the_solutions_inside_the_window(self):
        self._runPatchdock(minValue=7.5)

        params = _read(os.path.join(self.workDir, 'used_params.txt')).splitlines()
        # PatchDock counts atoms from 1, the state stores 0-based anchors.
        self.assertEqual(params[-1], 'distanceConstraints 16 25 14.0')
        self.assertIn('clusterParams 0.1 4 2.0 2.0', params)

        # Solution 1 (7.41 A) is below the minimum; 2 and 3 become pd.1 and pd.2.
        self.assertEqual(_read(os.path.join(self.workDir, 'transOutput_calls.txt')),
                         '2 3\n')
        results = os.path.join(self.workDir, 'Patchdock_Results')
        self.assertEqual(sorted(os.listdir(results)), ['pd.1.pdb', 'pd.2.pdb'])
        pd1 = _read(os.path.join(results, 'pd.1.pdb'))
        self.assertIn('solution 2', pd1)
        self.assertIn('PT1 Y', pd1)
        self.assertNotIn('PT1 X', pd1)
        self.assertIn('solution 3', _read(os.path.join(results, 'pd.2.pdb')))
        self.assertEqual(self._readState()['numResults'], 2)
        # Intermediate models are not left behind.
        self.assertEqual(glob.glob(os.path.join(self.workDir, 'Patchdock_output.*.pdb')), [])

    def test_runPatchdock_fails_when_no_solution_fits(self):
        with self.assertRaises(RuntimeError):
            self._runPatchdock(minValue=20.0, maxValue=30.0)

    # ------------------------------ other phases ------------------------------
    def _runWithStubs(self, func, args, **modules):
        with mock.patch.dict(sys.modules, modules), _chdir(self.workDir):
            return func(args)

    def test_runSampleDist_stores_the_window(self):
        self._writeState({'heads': ['Head0_H.sdf', 'Head1_H.sdf'], 'anchors': [3, 4]})
        calls = []

        def sampleDist(heads, anchors, smi):
            calls.append((heads, anchors, _read(smi)))
            return 7.5, 14.0

        self._runWithStubs(self.driver.runSampleDist, argparse.Namespace(smiles=' CCO '),
                           protac_lib=types.SimpleNamespace(SampleDist=sampleDist))
        self.assertEqual(calls, [(['Head0_H.sdf', 'Head1_H.sdf'], [3, 4], 'CCO\n')])
        state = self._readState()
        self.assertEqual((state['minValue'], state['maxValue']), (7.5, 14.0))

    def test_runSampleDist_fails_on_both_empty_results(self):
        for result in ((None, None), (0, 0)):
            self._writeState({'heads': ['a', 'b'], 'anchors': [0, 0]})
            stub = types.SimpleNamespace(SampleDist=lambda *a, r=result: r)
            with self.assertRaises(RuntimeError):
                self._runWithStubs(self.driver.runSampleDist,
                                   argparse.Namespace(smiles='CCO'), protac_lib=stub)

    def test_runLocalDocking_builds_one_command_per_solution(self):
        self._writeState({'ptParams': ['PT0.params', 'PT1.params'], 'numResults': 2})
        os.makedirs(os.path.join(self.workDir, 'Patchdock_Results'))
        _write(os.path.join(self.workDir, 'Patchdock_Results', 'local.fasc'), 'stale')
        captured = []
        rosettaStub = types.SimpleNamespace(
            local_docking=lambda pdb, c1, c2, p1, p2, n: ' '.join([pdb, c1, c2, p1, p2,
                                                                    str(n)]))
        args = argparse.Namespace(chain1='AC', chain2='B', nstruct=10, threads=3)
        with mock.patch.object(self.driver, '_runCommandsInParallel',
                               lambda cmds, cwd, threads, what: captured.append(
                                   (cmds, cwd, threads))):
            self._runWithStubs(self.driver.runLocalDocking, args, rosetta=rosettaStub)

        params0 = os.path.join(os.path.realpath(self.workDir), 'PT0.params')
        params1 = os.path.join(os.path.realpath(self.workDir), 'PT1.params')
        commands, cwd, threads = captured[0]
        self.assertEqual(commands, [f'pd.1.pdb ACX BY {params0} {params1} 10',
                                    f'pd.2.pdb ACX BY {params0} {params1} 10'])
        self.assertEqual((cwd, threads), ('Patchdock_Results', 3))
        self.assertFalse(os.path.exists(
            os.path.join(self.workDir, 'Patchdock_Results', 'local.fasc')))

    def test_runConstraintConf_builds_one_command_per_docking_model(self):
        self._writeState({'heads': ['Head0_H.sdf', 'Head1_H.sdf']})
        results = os.path.join(self.workDir, 'Patchdock_Results')
        for name in ('pd.1_docking_0001.pdb', 'pd.12_docking_0010.pdb', 'score.sc'):
            _write(os.path.join(results, name), 'x')
        captured = []
        args = argparse.Namespace(chain1='A', chain2='B', threads=2)
        with mock.patch.object(self.driver, '_runCommandsInParallel',
                               lambda cmds, cwd, threads, what: captured.append(cmds)):
            self._runWithStubs(self.driver.runConstraintConf, args)

        script = os.path.join(self._scriptsFol, 'constraint_generation.py')
        # Not a bare 'python': the Python 2 env is on PATH too.
        self.assertEqual(sorted(captured[0]), [
            f'{sys.executable} {script} ../Head0_H.sdf ../Head1_H.sdf ../protac.smi '
            f'12_10 pd.12_docking_0010.pdb AB',
            f'{sys.executable} {script} ../Head0_H.sdf ../Head1_H.sdf ../protac.smi '
            f'1_1 pd.1_docking_0001.pdb AB'])
        # score.sc is appended to by every job, so the old one must go.
        self.assertFalse(os.path.exists(os.path.join(results, 'score.sc')))

    def test_runConstraintConf_fails_without_docking_models(self):
        self._writeState({'heads': ['a', 'b']})
        os.makedirs(os.path.join(self.workDir, 'Patchdock_Results'))
        with self.assertRaises(RuntimeError):
            self._runWithStubs(self.driver.runConstraintConf,
                               argparse.Namespace(chain1='A', chain2='B', threads=1))

    def _runClustering(self, main):
        for name in ('Init0.pdb', 'Init1.pdb'):
            _write(os.path.join(self.workDir, name), name + '\n')
        os.makedirs(os.path.join(self.workDir, 'Patchdock_Results'))
        args = argparse.Namespace(chain2='B', topScore=1000, topLocal=200, rmsd=4.0)
        cwd = os.getcwd()
        try:
            self._runWithStubs(self.driver.runClustering, args,
                               clustering=types.SimpleNamespace(main=main))
        finally:
            # The driver must bring the cwd back even when clustering moved it.
            self.assertEqual(os.getcwd(), cwd)

    def test_runClustering_passes_its_arguments(self):
        calls = []

        def main(name, argv):
            calls.append(argv)
            os.chdir('Patchdock_Results')
            os.mkdir('../Results')

        self._runClustering(main)
        self.assertEqual(calls, [['1000', '200', '4.0', 'B']])
        self.assertEqual(_read(os.path.join(self.workDir, 'Init.pdb')),
                         'Init0.pdb\nInit1.pdb\n')

    def test_runClustering_turns_sys_exit_into_an_error(self):
        # clustering.py calls a bare sys.exit(), which would leave the step green.
        def main(name, argv):
            os.chdir('Patchdock_Results')
            sys.exit()

        with self.assertRaises(RuntimeError):
            self._runClustering(main)

    def test_runClustering_fails_without_results(self):
        # clustering.main() just returns when score.sc is missing or no model is < 0.
        def main(name, argv):
            os.chdir('Patchdock_Results')

        with self.assertRaises(RuntimeError):
            self._runClustering(main)

    # ------------------------------ parallel runner ------------------------------
    def test_runCommandsInParallel_runs_every_command(self):
        commands = [f'touch out{i}' for i in range(5)]
        self.driver._runCommandsInParallel(commands, self.workDir, 2, 'test jobs')
        for i in range(5):
            self.assertTrue(os.path.exists(os.path.join(self.workDir, f'out{i}')))

    def test_runCommandsInParallel_tolerates_some_failures(self):
        # Same as the original: each command is an independent cluster job, and the
        # later phases work with whatever models exist.
        self.driver._runCommandsInParallel(['touch survivor', 'exit 1'],
                                           self.workDir, 2, 'test jobs')
        self.assertTrue(os.path.exists(os.path.join(self.workDir, 'survivor')))

    def test_runCommandsInParallel_raises_when_all_fail(self):
        with self.assertRaises(RuntimeError):
            self.driver._runCommandsInParallel(['exit 1', 'exit 2'], self.workDir, 2,
                                               'test jobs')

    # ------------------------------ warhead normalisation ------------------------------
    @unittest.skipUnless(importlib.util.find_spec('rdkit'),
                         'RDKit is not installed in this environment; this test only runs '
                         'where RDKit is, e.g. the PRosettaC Python env')
    def test_normalizeHead_strips_H_and_remaps_the_anchor(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles('Cc1ccccc1O'))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        # Put a hydrogen before the O so the anchor index must shift once H are gone.
        firstH = next(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        order = [firstH] + [i for i in range(mol.GetNumAtoms()) if i != firstH]
        mol = Chem.RenumberAtoms(mol, order)
        sdf = os.path.join(self.workDir, 'head.sdf')
        Chem.MolToMolFile(mol, sdf, kekulize=True)
        oIdx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'O')

        newAnchor = self.driver._normalizeHead(sdf, oIdx)

        out = Chem.SDMolSupplier(sdf, removeHs=False, sanitize=False)[0]
        self.assertTrue(all(a.GetAtomicNum() != 1 for a in out.GetAtoms()))
        self.assertEqual(out.GetAtomWithIdx(newAnchor).GetSymbol(), 'O')
        with self.assertRaises(RuntimeError):
            self.driver._normalizeHead(sdf, out.GetNumAtoms())
        with self.assertRaises(RuntimeError):
            self.driver._normalizeHead(sdf, -1)


class TestPRosettaCProtocol(BaseTest):
    """ The protocol's own code: validation, input conversion and output parsing. Steps
    are called directly on a saved protocol, so no external tool is needed. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        protImport = cls.newProtocol(ProtImportPdb, inputPdbData=1, pdbFile=_TOY_PDB)
        cls.launchProtocol(protImport)
        cls.structure = protImport.outputPdb

        # Real SetOfSmallMolecules, attached as outputs of the import run so that the
        # protocol can point at them.
        dataDir = protImport._getExtraPath()
        cls.headFiles = {name: _write(os.path.join(dataDir, f'{name}.sdf'), f'{name}\n')
                         for name in ('lig1', 'lig2', 'lig2_pose', 'dup_a', 'dup_b')}
        heads = SetOfSmallMolecules().create(dataDir, suffix='heads')
        heads.append(SmallMolecule(smallMolFilename=cls.headFiles['lig1'], molName='lig1'))
        heads.append(SmallMolecule(smallMolFilename=cls.headFiles['lig2'], molName='lig2',
                                   poseFile=cls.headFiles['lig2_pose']))
        heads.append(SmallMolecule(smallMolFilename=cls.headFiles['lig1'], molName='lig3'))
        # Two poses of one molecule, as a docking run produces them.
        dupHeads = SetOfSmallMolecules().create(dataDir, suffix='dup')
        dupHeads.append(SmallMolecule(smallMolFilename=cls.headFiles['dup_a'],
                                      molName='dup', poseId=1))
        dupHeads.append(SmallMolecule(smallMolFilename=cls.headFiles['dup_b'],
                                      molName='dup', poseId=2))
        protImport._defineOutputs(heads=heads, dupHeads=dupHeads)
        protImport._store()
        cls.heads = protImport.heads
        cls.dupHeads = protImport.dupHeads

    def _newProt(self, **overrides):
        args = {'structure1': self.structure, 'structure2': self.structure,
                'heads1': self.heads, 'heads2': self.heads,
                'chain1': 'A', 'chain2': 'B', 'head1Name': 'lig1', 'head2Name': 'lig2',
                'anchor1': 1, 'anchor2': 1, 'protacSmiles': 'CCO'}
        args.update(overrides)
        return self.newProtocol(ProtPRosettaC, **args)

    def _newSavedProt(self, **overrides):
        prot = self._newProt(**overrides)
        self.saveProtocol(prot)
        prot.makeWorkingDir()
        return prot

    # ------------------------------ registration ------------------------------
    def test_protocol_is_registered(self):
        from protac import protocols
        self.assertIs(getattr(protocols, 'ProtPRosettaC', None), ProtPRosettaC)
        import protac
        with open(os.path.join(os.path.dirname(protac.__file__), 'protocols.conf')) as f:
            self.assertIn('ProtPRosettaC', f.read())

    # -------------------------------- validation -------------------------------
    def test_validate_accepts_a_correct_form(self):
        self.assertEqual(self._newProt(chain1='AC')._validate(), [])

    def _assertInvalid(self, fragment, **overrides):
        errors = self._newProt(**overrides)._validate()
        self.assertTrue(any(fragment in e for e in errors), errors)

    def test_validate_rejects_overlapping_chains(self):
        self._assertInvalid('must not overlap', chain1='AC', chain2='C')

    def test_validate_rejects_a_multi_chain_moving_structure(self):
        # clustering.py compares the PDB chain column character by character.
        self._assertInvalid('single chain', chain2='BC')

    def test_validate_rejects_lowercase_and_reserved_chains(self):
        # clean_pdb.py upper-cases the chains; X and Y are the warheads' chains.
        for chain1 in ('a', 'X', 'Y', 'A-'):
            self._assertInvalid('uppercase letters or digits', chain1=chain1)

    def test_validate_rejects_non_positive_sizes_and_anchors(self):
        for field in ('anchor1', 'anchor2', 'patchdockResults', 'localNstruct'):
            self._assertInvalid('at least 1', **{field: 0})

    def test_validate_rejects_more_local_than_total_models(self):
        self._assertInvalid('cannot be larger', clusterTopScore=100, clusterTopLocal=200)

    def test_validate_rejects_a_warhead_name_missing_from_its_set(self):
        self._assertInvalid('no molecule named "nope"', head2Name='nope')

    def test_validate_rejects_an_ambiguous_warhead_name(self):
        self._assertInvalid('2 molecules are named "dup"', heads1=self.dupHeads,
                            head1Name='dup')

    def test_validate_reports_a_missing_installation_once(self):
        prot = self._newProt()
        errors = prot.validate()
        # Protocol.validate() collects validateInstallation() by itself, so _validate()
        # must not repeat it.
        for installError in prot.validateInstallation():
            self.assertEqual(errors.count(installError), 1)

    # ------------------------------ step arguments ------------------------------
    def _stepArgs(self, **overrides):
        prot = self._newProt(**overrides)
        prot._insertAllSteps()
        return [(step.funcName.get(), step.argsStr.get()) for step in prot._steps]

    def test_changing_an_input_reruns_its_step_and_the_later_ones(self):
        # Scipion only reruns a step whose arguments changed when continuing a run.
        base = self._stepArgs()
        names = [name for name, _ in base]
        self.assertEqual(names, ['convertInputStep', 'prepareStructuresStep',
                                 'sampleDistStep', 'patchdockStep', 'localDockingStep',
                                 'constraintConfStep', 'clusteringStep', 'createOutputStep'])
        for field, value, firstChanged in (('head1Name', 'lig3', 0), ('anchor1', 2, 1),
                                           ('protacSmiles', 'CCN', 2),
                                           ('patchdockResults', 7, 3),
                                           ('localNstruct', 3, 4),
                                           ('clusterRmsd', 3.0, 6)):
            changed = self._stepArgs(**{field: value})
            for i, (old, new) in enumerate(zip(base, changed)):
                if i < firstChanged:
                    self.assertEqual(old, new, (field, old[0]))
                else:
                    self.assertNotEqual(old, new, (field, old[0]))

    # ------------------------------ warhead lookup ------------------------------
    def test_findMolByName_returns_the_named_molecule(self):
        found = ProtPRosettaC._findMolByName(self.heads, 'lig2')
        self.assertEqual(found.getMolName(), 'lig2')
        self.assertEqual(found.getPoseFile(), self.headFiles['lig2_pose'])

    def test_findMolByName_raises_when_absent(self):
        with self.assertRaises(ValueError):
            ProtPRosettaC._findMolByName(self.heads, 'nope')

    def test_findMolByName_rejects_an_ambiguous_name(self):
        # Poses of a docked set share their name: picking one silently would hide it.
        with self.assertRaises(ValueError):
            ProtPRosettaC._findMolByName(self.dupHeads, 'dup')

    def test_convertInputStep_prefers_the_pose_file(self):
        prot = self._newSavedProt()
        prot.convertInputStep('lig1', 'lig2', 0, 0)
        self.assertEqual(_read(prot._getHeadFile(1)), 'lig1\n')
        self.assertEqual(_read(prot._getHeadFile(2)), 'lig2_pose\n')

    # ------------------------------ score file ------------------------------
    def _writeScoreFile(self, content):
        return _write(os.path.join(tempfile.mkdtemp(prefix='prosettac_score_'), 'score.sc'),
                      content)

    def test_parseScoreFile_reads_columns_by_name(self):
        scoreFile = self._writeScoreFile(
            'SEQUENCE: \n'
            'SCORE: total_score fa_atr description\n'
            'SCORE: -12.5 -3.0 combined_1_1_0001\n'
            'SCORE: -8.25 -2.0 combined_2_3_0001\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile),
                         {'combined_1_1_0001': -12.5, 'combined_2_3_0001': -8.25})

    def test_parseScoreFile_survives_reordered_and_extra_columns(self):
        # A Rosetta build with extra score terms must not shift what gets read.
        scoreFile = self._writeScoreFile('SCORE: description extra_term total_score\n'
                                         'SCORE: combined_1_1_0001 0.1 -12.5\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile),
                         {'combined_1_1_0001': -12.5})

    def test_parseScoreFile_skips_unparsable_lines(self):
        scoreFile = self._writeScoreFile('SCORE: total_score description\n'
                                         'not a score line\n'
                                         'SCORE: notanumber model_bad\n'
                                         'SCORE: total_score description\n'
                                         'SCORE: -1.0 model_good\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile), {'model_good': -1.0})

    def test_parseScoreFile_returns_empty_when_missing(self):
        # A missing score file is metadata loss, not a reason to fail a finished run.
        missing = os.path.join(tempfile.mkdtemp(), 'score.sc')
        self.assertEqual(ProtPRosettaC._parseScoreFile(missing), {})

    # ------------------------------ output ------------------------------
    def test_getClusterDirs_sorts_numerically_and_ignores_other_entries(self):
        prot = self._newSavedProt()
        for name in ('cluster1', 'cluster2', 'cluster10', 'clusterX'):
            os.makedirs(os.path.join(prot._getResultsDir(), name))
        _write(os.path.join(prot._getResultsDir(), 'cluster3'), 'a file')
        # A plain sort would put cluster10 before cluster2.
        self.assertEqual([clusterId for clusterId, _ in prot._getClusterDirs()], [1, 2, 10])

    def test_createOutputStep_one_item_per_clustered_model(self):
        prot = self._newSavedProt()
        # Same names as PRosettaC: combined_<solution>_<model>_0001 in score.sc and Results/.
        _write(os.path.join(prot._getWorkDirFile('Patchdock_Results'), 'score.sc'),
               'SEQUENCE: \n'
               'SCORE: total_score description\n'
               'SCORE: -30.0 combined_1_1_0001\n'
               'SCORE: -20.0 combined_2_4_0001\n')
        pdb = _read(_TOY_PDB)
        _write(os.path.join(prot._getResultsDir(), 'cluster1', 'combined_1_1_0001.pdb'), pdb)
        _write(os.path.join(prot._getResultsDir(), 'cluster2', 'combined_2_4_0001.pdb'), pdb)
        _write(os.path.join(prot._getResultsDir(), 'cluster2', 'combined_3_1_0001.pdb'), pdb)

        prot.createOutputStep('[]')

        models = {os.path.basename(m.getFileName()): (m._clusterId.get(), m._score.get())
                  for m in prot.outputTernaryModels}
        self.assertEqual(models, {'combined_1_1_0001.pdb': (1, -30.0),
                                  'combined_2_4_0001.pdb': (2, -20.0),
                                  'combined_3_1_0001.pdb': (2, None)})

    def test_createOutputStep_fails_without_clusters(self):
        prot = self._newSavedProt()
        os.makedirs(prot._getResultsDir())
        with self.assertRaises(RuntimeError):
            prot.createOutputStep('[]')

    def test_createOutputStep_fails_on_empty_clusters(self):
        prot = self._newSavedProt()
        os.makedirs(os.path.join(prot._getResultsDir(), 'cluster1'))
        with self.assertRaises(RuntimeError):
            prot.createOutputStep('[]')
