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

Everything here runs without PatchDock/Rosetta/PRosettaC installed, on purpose:
scipion_testrunner only skips tests on missing *Python modules*, never on missing
binaries or *_HOME variables, so an environment-dependent test has to skip itself.
The single test that does need a real installation
(TestPRosettaCProtocol.test_installation_is_usable) does exactly that.

What is NOT covered here, and why: an end-to-end run of the protocol. It would need a
real pair of structures with their bound warheads plus a matching PROTAC SMILES, and
hours of PatchDock + Rosetta - unlike ProtPROTACTransplantWarhead, there is no synthetic
fixture that produces a meaningful result.
"""

import glob
import importlib.util
import os
import shutil
import sys
import tempfile
import unittest

from pyworkflow.tests import BaseTest, setupTestProject

from protac.protocols import ProtPRosettaC

_SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'scripts', 'run_prosettac.py')


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


class TestPRosettaCDriver(unittest.TestCase):
    """ The driver's own helpers - the part that is ours rather than PRosettaC's, and the
    only part that can be tested without the external tools. """

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
        # Both the sys.path entry the driver inserts at import time and the SCRIPTS_FOL
        # that _loadDriver sets for it point into the temporary tree removed below;
        # leaving either behind would follow this process into every later test.
        sys.path[:] = cls._pathBefore
        if cls._scriptsFolBefore is None:
            os.environ.pop('SCRIPTS_FOL', None)
        else:
            os.environ['SCRIPTS_FOL'] = cls._scriptsFolBefore
        shutil.rmtree(cls._tmpRoot, ignore_errors=True)

    def setUp(self):
        self.workDir = tempfile.mkdtemp(prefix='prosettac_case_', dir=self._tmpRoot)

    def test_dockingSuffix_matches_the_originals_own_format(self):
        # main.py derives this with "s.split('.')[1].split('_')"; clustering.py parses the
        # two numbers back out, so the leading zeros of the docking index must be dropped.
        self.assertEqual(self.driver._dockingSuffix('pd.7_docking_0003.pdb'), '7_3')
        self.assertEqual(self.driver._dockingSuffix('pd.12_docking_0050.pdb'), '12_50')

    def test_stageSimpleName_copies_under_a_dot_free_name(self):
        src = os.path.join(self.workDir, 'a.b.c', 'input.struct.pdb')
        os.makedirs(os.path.dirname(src))
        with open(src, 'w') as f:
            f.write('ATOM\n')

        cwd = os.getcwd()
        os.chdir(self.workDir)
        try:
            self.assertEqual(self.driver._stageSimpleName(src, 'Struct0.pdb'),
                             'Struct0.pdb')
            self.assertTrue(os.path.exists('Struct0.pdb'))
            with open('Struct0.pdb') as f:
                self.assertEqual(f.read(), 'ATOM\n')
        finally:
            os.chdir(cwd)

    def test_catFiles_concatenates_in_order(self):
        parts = []
        for i, content in enumerate(('first\n', 'second\n')):
            path = os.path.join(self.workDir, f'part{i}.pdb')
            with open(path, 'w') as f:
                f.write(content)
            parts.append(path)

        dest = os.path.join(self.workDir, 'joined.pdb')
        self.driver._catFiles(parts, dest)
        with open(dest) as f:
            # Order matters: it fixes the atom numbering PatchDock later sees.
            self.assertEqual(f.read(), 'first\nsecond\n')

    def test_removeStale_only_removes_matching_files(self):
        for name in ('local.fasc', 'pd.1_docking_0001.pdb', 'keepme.pdb'):
            with open(os.path.join(self.workDir, name), 'w') as f:
                f.write('')

        self.driver._removeStale(self.workDir, 'local.fasc', '*_docking_????.pdb')
        remaining = sorted(os.path.basename(p)
                           for p in glob.glob(os.path.join(self.workDir, '*')))
        self.assertEqual(remaining, ['keepme.pdb'])

    def test_removeStale_tolerates_a_missing_directory(self):
        # Phase 5/6 clear leftovers before checking anything exists.
        self.driver._removeStale(os.path.join(self.workDir, 'nope'), '*.pdb')

    def test_runCommandsInParallel_runs_every_command(self):
        commands = [f'touch out{i}' for i in range(5)]
        self.driver._runCommandsInParallel(commands, self.workDir, 2, 'test jobs')
        for i in range(5):
            self.assertTrue(os.path.exists(os.path.join(self.workDir, f'out{i}')))

    def test_runCommandsInParallel_tolerates_some_failures(self):
        # Matches the original's semantics: each command is an independent cluster job,
        # and the later phases work with whatever models exist.
        self.driver._runCommandsInParallel(['touch survivor', 'exit 1'],
                                           self.workDir, 2, 'test jobs')
        self.assertTrue(os.path.exists(os.path.join(self.workDir, 'survivor')))

    def test_runCommandsInParallel_raises_when_all_fail(self):
        with self.assertRaises(RuntimeError):
            self.driver._runCommandsInParallel(['exit 1', 'exit 2'], self.workDir, 2,
                                               'test jobs')

    @unittest.skipUnless(importlib.util.find_spec('rdkit'), 'needs RDKit (driver env)')
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

    def test_parseArgs_accepts_every_phase(self):
        argvBefore = list(sys.argv)
        try:
            for phase in ('prepare', 'sampledist', 'patchdock', 'localdocking',
                          'constraintconf', 'clustering'):
                sys.argv = ['run_prosettac.py', '--phase', phase]
                self.assertEqual(self.driver.parseArgs().phase, phase)
        finally:
            sys.argv = argvBefore


class _FakeMol:
    """ Minimal stand-in for pwchem's SmallMolecule: _findMolByName only calls
    getMolName()/clone() on it. """
    def __init__(self, name):
        self._name = name

    def getMolName(self):
        return self._name

    def clone(self):
        return self


class TestPRosettaCProtocol(BaseTest):
    """ The protocol's own pure helpers, its validation and its registration. Nothing here
    launches the protocol, so no external tool is needed. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @staticmethod
    def _newFormArgs(**overrides):
        args = {'chain1': 'A', 'chain2': 'B', 'head1Name': 'lig1', 'head2Name': 'lig2',
                'anchor1': 1, 'anchor2': 1, 'protacSmiles': 'CCO'}
        args.update(overrides)
        return args

    # ------------------------------ registration ------------------------------
    def test_protocol_is_registered(self):
        from protac import protocols
        self.assertIs(getattr(protocols, 'ProtPRosettaC', None), ProtPRosettaC)
        self.assertTrue(ProtPRosettaC._label)

    def test_protocol_is_in_the_scipion_menu(self):
        import protac
        confFile = os.path.join(os.path.dirname(protac.__file__), 'protocols.conf')
        with open(confFile) as f:
            self.assertIn('ProtPRosettaC', f.read())

    # -------------------------------- _validate -------------------------------
    def test_validate_rejects_overlapping_chains(self):
        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs(chain1='AC',
                                                                  chain2='C'))
        self.assertTrue(any('must not overlap' in e for e in prot._validate()))

    def test_validate_accepts_disjoint_chains(self):
        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs(chain1='AC',
                                                                  chain2='B'))
        self.assertFalse(any('must not overlap' in e for e in prot._validate()))

    def test_validate_rejects_a_multi_chain_moving_structure(self):
        # clustering.py compares the PDB chain column character by character, so this
        # would only blow up after every docking phase had already run.
        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs(chain1='A',
                                                                  chain2='BC'))
        self.assertTrue(any('single chain' in e for e in prot._validate()))

    def test_validate_reports_a_missing_installation(self):
        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs())
        try:
            prot._requireEnv()
        except FileNotFoundError as e:
            # Reported through validateInstallation(), which Protocol.validate() collects
            # by itself; _validate() must not repeat it or the form shows it twice.
            self.assertIn(str(e), prot.validateInstallation())
            self.assertNotIn(str(e), prot._validate())
        else:
            self.assertEqual(prot.validateInstallation(), [])

    def test_installation_is_usable(self):
        """ The only environment-dependent test: skips itself when the tools are not
        installed, since the test runner cannot skip on a *_HOME variable. """
        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs())
        errors = prot.validateInstallation()
        if errors:
            self.skipTest(f'PRosettaC/PatchDock/Rosetta not installed here: {errors[0]}')
        self.assertEqual(errors, [])

    # ------------------------------ pure helpers ------------------------------
    def test_findMolByName_returns_the_named_molecule(self):
        molecules = [_FakeMol('lig1'), _FakeMol('lig2')]
        found = ProtPRosettaC._findMolByName(molecules, 'lig2')
        self.assertEqual(found.getMolName(), 'lig2')

    def test_findMolByName_raises_when_absent(self):
        with self.assertRaises(ValueError):
            ProtPRosettaC._findMolByName([_FakeMol('lig1')], 'nope')

    def test_parseScoreFile_reads_columns_by_name(self):
        scoreFile = self._writeScoreFile(
            'SCORE: total_score fa_atr description\n'
            'SCORE: -12.5 -3.0 model_1\n'
            'SCORE: -8.25 -2.0 model_2\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile),
                         {'model_1': -12.5, 'model_2': -8.25})

    def test_parseScoreFile_survives_reordered_and_extra_columns(self):
        # The point of looking columns up by name: a Rosetta build with extra score terms
        # must not shift what gets read.
        scoreFile = self._writeScoreFile(
            'SCORE: description extra_term total_score\n'
            'SCORE: model_1 0.1 -12.5\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile), {'model_1': -12.5})

    def test_parseScoreFile_skips_unparsable_lines(self):
        scoreFile = self._writeScoreFile(
            'SCORE: total_score description\n'
            'not a score line\n'
            'SCORE: notanumber model_bad\n'
            'SCORE: -1.0 model_good\n')
        self.assertEqual(ProtPRosettaC._parseScoreFile(scoreFile), {'model_good': -1.0})

    def test_parseScoreFile_returns_empty_when_missing(self):
        # A missing score file is metadata loss, not a reason to fail a finished run.
        missing = os.path.join(tempfile.mkdtemp(), 'score.sc')
        self.assertEqual(ProtPRosettaC._parseScoreFile(missing), {})

    def test_getClusterDirs_sorts_numerically_and_ignores_other_entries(self):
        resultsDir = tempfile.mkdtemp(prefix='prosettac_results_')
        for name in ('cluster1', 'cluster2', 'cluster10', 'clusterX'):
            os.makedirs(os.path.join(resultsDir, name))
        with open(os.path.join(resultsDir, 'cluster3'), 'w') as f:
            f.write('a file, not a cluster directory')

        prot = self.newProtocol(ProtPRosettaC, **self._newFormArgs())
        prot._getResultsDir = lambda: resultsDir

        # A plain sort would put cluster10 before cluster2, silently mislabelling the
        # cluster ranking on every output model.
        self.assertEqual([clusterId for clusterId, _ in prot._getClusterDirs()],
                         [1, 2, 10])

    def _writeScoreFile(self, content):
        scoreFile = os.path.join(tempfile.mkdtemp(prefix='prosettac_score_'), 'score.sc')
        with open(scoreFile, 'w') as f:
            f.write(content)
        return scoreFile
