# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Name:     test of the plugin helpers in protac/__init__.py
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
Tests for the shims and environment helpers of the plugin class. The external tools are
replaced by fake directory trees and small shell scripts, so nothing needs installing.
"""

import os
import shutil
import stat
import subprocess
import tempfile
import unittest
from unittest import mock

from pwchem import Plugin as pwchemPlugin
from rosetta import Plugin as RosettaPlugin

from protac import Plugin
from protac.constants import (FRODOCK_BINARIES, FRODOCK_DIC, VINA_DIC, ADFRSUITE_DIC,
                              VOROMQA_DIC, FCC_DIC, PROTAC_MODEL_DIC,
                              PROTAC_MODEL_PYTHON_DIC, PROSETTAC_PYTHON2_DIC)
from protac.protocols import ProtPROTACModel
from protac.utils.molecules import ligandCentroid, listMolNames


def _touch(path, content='', executable=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    if executable:
        os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


class TestPluginHelpers(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='protac_plugin_')
        self.addCleanup(shutil.rmtree, self.tmp, True)

    # ------------------------------ tool homes ------------------------------
    def test_requireToolHome(self):
        home = os.path.join(self.tmp, 'vina')
        _touch(os.path.join(home, 'bin', 'vina'))
        for value in (None, os.path.join(self.tmp, 'missing')):
            with mock.patch.object(Plugin, 'getVar', return_value=value), \
                    self.assertRaises(FileNotFoundError):
                Plugin._requireToolHome(VINA_DIC)
        with mock.patch.object(Plugin, 'getVar', return_value=home):
            self.assertEqual(Plugin._requireToolHome(VINA_DIC, 'bin/vina'), home)
            with self.assertRaises(FileNotFoundError):
                Plugin._requireToolHome(VINA_DIC, 'bin/nope')

    def test_requirePatchdockHome(self):
        for value in (None, os.path.join(self.tmp, 'missing')):
            with mock.patch.object(Plugin, 'getVar', return_value=value), \
                    self.assertRaises(FileNotFoundError):
                Plugin._requirePatchdockHome()
        with mock.patch.object(Plugin, 'getVar', return_value=self.tmp):
            self.assertEqual(Plugin._requirePatchdockHome(), self.tmp)

    # ------------------------------ FRODOCK ------------------------------
    def _frodockHome(self):
        home = os.path.join(self.tmp, 'frodock')
        for name in FRODOCK_BINARIES:
            _touch(os.path.join(home, 'bin', name))
            _touch(os.path.join(home, 'bin', f'{name}_gcc'))
        _touch(os.path.join(home, 'bin', 'soap.bin'))
        return home

    def _prepareFrodock(self, home, loadable):
        shim = os.path.join(self.tmp, 'shim')
        with mock.patch.object(Plugin, '_requireToolHome', return_value=home), \
                mock.patch.object(Plugin, '_binaryLoads', side_effect=lambda p: p in loadable):
            Plugin.prepareFrodockBinDir(shim)
        return shim

    def test_prepareFrodockBinDir_prefers_intel_and_falls_back_to_gcc(self):
        home = self._frodockHome()
        binDir = os.path.join(home, 'bin')
        # Intel loads for every program except frodock itself.
        loadable = {os.path.join(binDir, n) for n in FRODOCK_BINARIES if n != 'frodock'}
        loadable.add(os.path.join(binDir, 'frodock_gcc'))
        shim = self._prepareFrodock(home, loadable)

        for name in FRODOCK_BINARIES:
            expected = f'{name}_gcc' if name == 'frodock' else name
            self.assertEqual(os.path.realpath(os.path.join(shim, 'bin', name)),
                             os.path.realpath(os.path.join(binDir, expected)))
        self.assertTrue(os.path.exists(os.path.join(shim, 'bin', 'soap.bin')))

    def test_prepareFrodockBinDir_fails_when_no_build_loads(self):
        with self.assertRaises(FileNotFoundError):
            self._prepareFrodock(self._frodockHome(), loadable=set())

    def test_prepareFrodockBinDir_needs_soap_bin(self):
        home = self._frodockHome()
        os.remove(os.path.join(home, 'bin', 'soap.bin'))
        loadable = {os.path.join(home, 'bin', n) for n in FRODOCK_BINARIES}
        with self.assertRaises(FileNotFoundError):
            self._prepareFrodock(home, loadable)

    # ------------------------------ Rosetta ------------------------------
    def _rosettaHome(self, builds, toolsUnderMain=False):
        home = os.path.join(self.tmp, 'rosetta')
        for build in builds:
            _touch(os.path.join(home, 'main', 'source', 'bin',
                                f'rosetta_scripts.{build}.linuxgccrelease'))
        _touch(os.path.join(home, 'main', 'source', 'scripts', 'python', 'public',
                            'molfile_to_params.py'))
        _touch(os.path.join(home, 'main', 'database', 'chemical', 'x'))
        tools = os.path.join(home, 'main', 'tools') if toolsUnderMain else os.path.join(
            home, 'tools')
        _touch(os.path.join(tools, 'protein_tools', 'scripts', 'clean_pdb.py'))
        return home

    def test_requireRosettaScriptsBinary_build_order(self):
        cases = [(['static', 'default', 'mpi'], 'default'), (['mpi', 'static'], 'static'),
                 (['mpi'], 'mpi'), (['cxx11thread', 'boost'], 'boost')]
        for builds, expected in cases:
            home = self._rosettaHome(builds)
            with mock.patch.object(RosettaPlugin, 'getVar', return_value=home):
                self.assertTrue(Plugin.requireRosettaScriptsBinary().endswith(
                    f'rosetta_scripts.{expected}.linuxgccrelease'), builds)
            shutil.rmtree(home)

    def test_requireRosettaScriptsBinary_fails_without_binary_or_home(self):
        home = self._rosettaHome([])
        for value in (None, os.path.join(self.tmp, 'missing'), home):
            with mock.patch.object(RosettaPlugin, 'getVar', return_value=value), \
                    self.assertRaises(FileNotFoundError):
                Plugin.requireRosettaScriptsBinary()

    def test_prepareRosettaScriptsShim_exposes_the_layout_prosettac_expects(self):
        # Binary bundles keep tools/ under main/, source checkouts at the top level.
        for toolsUnderMain in (False, True):
            home = self._rosettaHome(['static'], toolsUnderMain=toolsUnderMain)
            shim = os.path.join(self.tmp, 'shim')
            with mock.patch.object(RosettaPlugin, 'getVar', return_value=home):
                Plugin.prepareRosettaScriptsShim(shim)
            for relPath in ('main/source/bin/rosetta_scripts.default.linuxgccrelease',
                            'main/source/scripts/python/public/molfile_to_params.py',
                            'tools/protein_tools/scripts/clean_pdb.py',
                            'main/database/chemical/x'):
                self.assertTrue(os.path.exists(os.path.join(shim, relPath)),
                                (toolsUnderMain, relPath))
            shutil.rmtree(home)
            shutil.rmtree(shim)

    def test_prepareRosettaScriptsShim_needs_the_database(self):
        home = self._rosettaHome(['static'])
        shutil.rmtree(os.path.join(home, 'main', 'database'))
        with mock.patch.object(RosettaPlugin, 'getVar', return_value=home), \
                self.assertRaises(FileNotFoundError):
            Plugin.prepareRosettaScriptsShim(os.path.join(self.tmp, 'shim'))

    # ------------------------------ OpenBabel ------------------------------
    def _babel(self, *args):
        """ Runs the generated babel wrapper against an obabel stand-in that logs its
        arguments and writes a marker plus its input to the -O file. """
        fakeObabel = _touch(os.path.join(self.tmp, 'obabel'),
                            '#!/bin/bash\n'
                            f'echo "$@" >> "{self.tmp}/obabel.log"\n'
                            '{ echo converted; cat "$1"; } > "$3"\n', executable=True)
        shim = os.path.join(self.tmp, 'babel_shim')
        with mock.patch.object(Plugin, 'requireObabelBinary', return_value=fakeObabel), \
                mock.patch.object(pwchemPlugin, 'getEnvActivationCommand',
                                  return_value='true'):
            Plugin.preparePRosettaCBabelShim(shim)
        subprocess.run([os.path.join(shim, 'babel')] + list(args), cwd=self.tmp,
                       check=True)
        with open(os.path.join(self.tmp, 'obabel.log')) as f:
            return f.read().splitlines()

    def test_babel_shim_translates_the_call(self):
        _touch(os.path.join(self.tmp, 'in.pdb'), 'ATOM\n')
        self.assertEqual(self._babel('in.pdb', 'out.sdf', '-h'), ['in.pdb -O out.sdf -h'])
        with open(os.path.join(self.tmp, 'out.sdf')) as f:
            self.assertEqual(f.read(), 'converted\nATOM\n')

    def test_babel_shim_handles_in_place_conversion(self):
        # PRosettaC converts some files onto themselves.
        _touch(os.path.join(self.tmp, 'head.sdf'), 'MOL\n')
        log = self._babel('head.sdf', 'head.sdf', '-h')
        self.assertEqual(len(log), 1)
        self.assertNotIn('-O head.sdf', log[0])
        with open(os.path.join(self.tmp, 'head.sdf')) as f:
            self.assertEqual(f.read(), 'converted\nMOL\n')

    # ------------------------------ environments ------------------------------
    def test_getPRosettaCEnviron(self):
        with mock.patch.object(Plugin, '_requireToolHome', return_value='/opt/py2'), \
                mock.patch.object(Plugin, '_requirePatchdockHome', return_value='/opt/pd'), \
                mock.patch.object(Plugin, 'getPRosettaCScript',
                                  return_value='/opt/prosettac/prosettac'), \
                mock.patch.dict(os.environ, {'PATH': '/usr/bin:/bin'}):
            env = Plugin.getPRosettaCEnviron('rosetta_shim', 'babel_shim')
            Plugin._requireToolHome.assert_called_with(PROSETTAC_PYTHON2_DIC)

        self.assertEqual(env['PATCHDOCK'], '/opt/pd')
        # rosetta.py appends file names to SCRIPTS_FOL without a separator.
        self.assertEqual(env['SCRIPTS_FOL'], '/opt/prosettac/prosettac/')
        self.assertEqual(env['ROSETTA_FOL'], os.path.abspath('rosetta_shim'))
        self.assertEqual(env['OB'], os.path.abspath('babel_shim'))
        # Python 2 goes last, or it would shadow the activated env's python.
        self.assertEqual(env['PATH'], '/usr/bin:/bin' + os.pathsep + '/opt/py2/bin')

    def _protacModelEnviron(self, rosettaHome='/opt/rosetta', **kwargs):
        # Rosetta's home must exist on disk to be accepted.
        with mock.patch('os.path.isdir', return_value=True), mock.patch.object(Plugin, '_requireToolHome',
                               side_effect=lambda dic, binary='': f"/opt/{dic['name']}"), \
                mock.patch.object(RosettaPlugin, 'getVar', return_value=rosettaHome):
            return Plugin.getProtacModelEnviron(**kwargs)

    def test_getProtacModelEnviron(self):
        expected = {'FRODOCK': '/opt/shim', 'ADFRSUITE': '/opt/adfrsuite',
                    'VINA': '/opt/vina', 'VOROMQA': '/opt/voromqa', 'FCC': '/opt/fcc/fcc',
                    'PROTAC_MODEL_HOME': '/opt/protac-model/protac-model'}
        self.assertEqual(self._protacModelEnviron(frodockHome='/opt/shim'), expected)
        self.assertEqual(self._protacModelEnviron(frodockHome='/opt/shim', needRosetta=True),
                         dict(expected, ROSETTA='/opt/rosetta'))
        self.assertEqual(self._protacModelEnviron()['FRODOCK'], f"/opt/{FRODOCK_DIC['name']}")

    def test_getProtacModelEnviron_needs_rosetta_only_to_refine(self):
        self.assertNotIn('ROSETTA', self._protacModelEnviron(rosettaHome=None))
        with self.assertRaises(FileNotFoundError):
            self._protacModelEnviron(rosettaHome=None, needRosetta=True)

    def test_getProtacModelEnviron_rejects_a_missing_rosetta_dir(self):
        missing = os.path.join(tempfile.gettempdir(), 'no_such_rosetta_dir')
        with mock.patch.object(Plugin, '_requireToolHome',
                               side_effect=lambda dic, binary='': f"/opt/{dic['name']}"), \
                mock.patch.object(RosettaPlugin, 'getVar', return_value=missing):
            with self.assertRaises(FileNotFoundError):
                Plugin.getProtacModelEnviron(needRosetta=True)

    # ------------------------------ PROTAC-Model installation ------------------------------
    def _protacModelInstall(self):
        """ A complete fake install, and the getVar() that points at it. """
        homes = {dic['home']: os.path.join(self.tmp, dic['name']) for dic in
                 (FRODOCK_DIC, ADFRSUITE_DIC, VINA_DIC, VOROMQA_DIC, FCC_DIC,
                  PROTAC_MODEL_DIC, PROTAC_MODEL_PYTHON_DIC)}
        files = [os.path.join(homes['ADFRSUITE_HOME'], 'bin', name) for name in
                 ('obabel', 'obenergy', 'prepare_ligand', 'prepare_receptor', 'reduce')]
        files += [os.path.join(homes['PROTAC_VINA_HOME'], 'bin', 'vina'),
                  os.path.join(homes['VOROMQA_HOME'], 'bin', 'voronota-voromqa'),
                  os.path.join(homes['FCC_HOME'], 'fcc', 'src', 'contact_fcc'),
                  os.path.join(homes['FCC_HOME'], 'fcc', 'scripts', 'make_contacts.py'),
                  os.path.join(homes['FRODOCK_HOME'], 'bin', 'soap.bin'),
                  os.path.join(homes['PROTAC_MODEL_HOME'], 'protac-model', 'utils',
                               'frodock.py'),
                  os.path.join(homes['PROTAC_MODEL_PYTHON_HOME'], 'bin', 'python2')]
        # One program only in its gcc build: either build is enough.
        files += [os.path.join(homes['FRODOCK_HOME'], 'bin',
                               name + ('_gcc' if name == 'frodock' else ''))
                  for name in FRODOCK_BINARIES]
        for path in files:
            _touch(path)
        return homes, mock.patch.object(Plugin, 'getVar', side_effect=homes.get)

    def test_checkProtacModelBinaries(self):
        homes, getVar = self._protacModelInstall()
        with getVar:
            self.assertEqual(Plugin.checkProtacModelBinaries(), [])
            self.assertEqual(ProtPROTACModel.validateInstallation(), [])

            os.remove(os.path.join(homes['PROTAC_VINA_HOME'], 'bin', 'vina'))
            os.remove(os.path.join(homes['FRODOCK_HOME'], 'bin', 'frodock_gcc'))
            errors = Plugin.checkProtacModelBinaries()
            self.assertEqual(len(errors), 2, errors)
            self.assertTrue(any('frodock_gcc' in e for e in errors))

            # A missing vina already fails the home check, with its own message.
            installErrors = ProtPROTACModel.validateInstallation()
            self.assertEqual(len(installErrors), 1)
            self.assertIn('PROTAC_VINA_HOME', installErrors[0])

    def test_validateInstallation_lists_missing_binaries(self):
        homes, getVar = self._protacModelInstall()
        os.remove(os.path.join(homes['ADFRSUITE_HOME'], 'bin', 'reduce'))
        with getVar:
            errors = ProtPROTACModel.validateInstallation()
        self.assertEqual(len(errors), 1)
        self.assertIn(os.path.join('bin', 'reduce'), errors[0])

    # ------------------------------ launching ------------------------------
    def test_runProgram_goes_through_the_protocol(self):
        calls = []
        protocol = mock.Mock(runJob=lambda *a, **kw: calls.append((a, kw)))
        with mock.patch.dict(os.environ, {'KEEP': '1', 'OVERRIDE': 'old'}):
            Plugin.runProgram(protocol, 'prog', '--x 1', extraEnvDict={'OVERRIDE': 'new'},
                              cwd='relative/dir')
        (args, kwargs), = calls
        self.assertEqual(args, ('prog', '--x 1'))
        self.assertEqual((kwargs['env']['KEEP'], kwargs['env']['OVERRIDE']), ('1', 'new'))
        self.assertEqual(kwargs['cwd'], os.path.abspath('relative/dir'))
        # An mpirun prefix would break the conda activation chain.
        self.assertEqual(kwargs['numberOfMpi'], 1)

    def test_runCondaScript_activates_the_env_first(self):
        protocol = mock.Mock()
        with mock.patch.object(Plugin, 'getEnvActivationCommand', return_value='activate x'):
            Plugin.runCondaScript(protocol, 'scripts/run.py', '--phase a', {'name': 'x'})
        program = protocol.runJob.call_args[0][0]
        self.assertEqual(program, f'activate x && python "{os.path.abspath("scripts/run.py")}"')


class _Mol:
    def __init__(self, name):
        self._name = name

    def getMolName(self):
        return self._name


class TestMoleculeUtils(unittest.TestCase):
    """ Helpers behind the wizards and the molecule lookups. """

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='protac_mols_')
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def _pdb(self, lines):
        return _touch(os.path.join(self.tmp, 'struct.pdb'), ''.join(
            '%-6s%5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n'
            % (rec, i, atom, res, chain, resi, x, y, z)
            for i, (rec, atom, res, chain, resi, x, y, z) in enumerate(lines, start=1)))

    def test_listMolNames_is_sorted_and_distinct(self):
        self.assertEqual(listMolNames([_Mol('b'), _Mol('a'), _Mol('b')]), ['a', 'b'])
        self.assertEqual(listMolNames(None), [])

    def test_ligandCentroid_averages_the_named_residue_only(self):
        pdb = self._pdb([('ATOM', 'CA', 'ALA', 'A', 1, 9.0, 9.0, 9.0),
                         ('HETATM', 'C1', 'LIG', 'A', 101, 0.0, 0.0, 0.0),
                         ('HETATM', 'C2', 'LIG', 'A', 101, 2.0, 4.0, -6.0),
                         ('HETATM', 'O', 'HOH', 'A', 201, 50.0, 50.0, 50.0)])
        self.assertEqual(ligandCentroid(pdb, 'lig'), (1.0, 2.0, -3.0))

    def test_ligandCentroid_rejects_missing_and_repeated_ligands(self):
        pdb = self._pdb([('HETATM', 'C1', 'LIG', 'A', 101, 0.0, 0.0, 0.0),
                         ('HETATM', 'C1', 'LIG', 'B', 101, 4.0, 0.0, 0.0)])
        with self.assertRaisesRegex(ValueError, 'no heteroatom'):
            ligandCentroid(pdb, 'XYZ')
        with self.assertRaisesRegex(ValueError, 'chains A, B'):
            ligandCentroid(pdb, 'LIG')

    def test_wizards_are_registered(self):
        from pwchem.wizards import SelectElementWizard, SelectLigandWizard
        from protac.protocols import ProtPRosettaC
        from protac.wizards import LigandCentroidWizard

        def targets(wizard, protocol):
            return {t for prot, ts in wizard._targets if prot is protocol for t in ts}

        self.assertEqual(targets(LigandCentroidWizard, ProtPROTACModel), {'siteCoords'})
        self.assertEqual(targets(SelectLigandWizard, ProtPROTACModel),
                         {'receptorLigandName', 'targetLigandName'})
        self.assertEqual(targets(SelectElementWizard, ProtPROTACModel),
                         {'e3Ligand1Name', 'e3Ligand2Name'})
        self.assertEqual(targets(SelectElementWizard, ProtPRosettaC),
                         {'head1Name', 'head2Name'})
