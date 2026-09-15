# -*- coding: utf-8 -*-
# **************************************************************************
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

import os
import subprocess

import pyworkflow.utils as pwutils
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from rosetta import Plugin as RosettaPlugin, ROSETTA_DIC

from .constants import *

# Computed once at import, before any chdir() the protocol's steps do later,
# abspath(__file__) could resolve wrong if computed after the cwd has moved.
_PLUGIN_DIR = os.path.dirname(os.path.abspath(__file__))

_version_ = "0.1"
# FRODOCK is a separate external tool, used by the PROTAC-Model pipeline
# for the initial global protein-protein docking step.
FRODOCK_DIC = {'name': 'frodock', 'home': 'FRODOCK_HOME'}

# Used later by filterPosesStep (PROTAC-Model's filter_frodock()): ADFRsuite binaries,
# Vina, Voromqa, and FCC's clustering scripts. No license click-through needed (unlike
# FRODOCK/Rosetta below), so these get real defineBinaries() via InstallHelper.
ADFRSUITE_DIC = {'name': 'adfrsuite', 'version': '1.0', 'home': 'ADFRSUITE_HOME'}
# Pinned for reproducibility, not as a workaround. The actual Vina bug (uninitialised
# grid box in --score_only, changed "Affinity" wording) is patched at runtime in
# run_protac_model.py, not here.
VINA_DIC = {'name': 'vina', 'version': '1.2.2', 'home': 'VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'version': '1.29.4816', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'version': 'latest', 'home': 'FCC_HOME'}

# PROTAC-Model's own code (main.py, utils/*) is called directly, not reimplemented
# Its code is genuine Python 2, so PROTAC_MODEL_PYTHON_HOME is a dedicated 
# Python 2.7+rdkit conda env, never scipion3's own.
PROTAC_MODEL_DIC = {'name': 'protac-model', 'version': 'latest', 'home': 'PROTAC_MODEL_HOME'}
PROTAC_MODEL_PYTHON_DIC = {'name': 'protac-model-python', 'version': '2.7',
                           'home': 'PROTAC_MODEL_PYTHON_HOME'}

class Plugin(pwchemPlugin):
    _homeVar = PROTAC_MODEL_DIC['home']
    _pathVars = [PROTAC_MODEL_DIC['home']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. """
        # FRODOCK_HOME stays manual.The user points it at their own download.
        cls._defineVar(FRODOCK_DIC['home'], None)

        # The five tools below have no license gate, _defineEmVar wires each home
        # to wherever defineBinaries()/InstallHelper installs it.
        cls._defineEmVar(ADFRSUITE_DIC['home'], cls.getEnvName(ADFRSUITE_DIC))
        cls._defineEmVar(VINA_DIC['home'], cls.getEnvName(VINA_DIC))
        cls._defineEmVar(VOROMQA_DIC['home'], cls.getEnvName(VOROMQA_DIC))
        cls._defineEmVar(FCC_DIC['home'], cls.getEnvName(FCC_DIC))
        cls._defineEmVar(PROTAC_MODEL_DIC['home'], cls.getEnvName(PROTAC_MODEL_DIC))
        cls._defineEmVar(PROTAC_MODEL_PYTHON_DIC['home'], cls.getEnvName(PROTAC_MODEL_PYTHON_DIC))

    @classmethod
    def defineBinaries(cls, env):
        # FRODOCK/Rosetta excluded on purpose (manual,ROSETTA_HOME 
        # is scipion-chem-rosetta's own responsibility).
        cls.addADFRSuitePackage(env)
        cls.addVinaPackage(env)
        cls.addVoromqaPackage(env)
        cls.addFCCPackage(env)
        cls.addProtacModelPackage(env)
        cls.addProtacModelPythonPackage(env)

    # ---------------------------- Package installers (InstallHelper) -------------
    @classmethod
    def addADFRSuitePackage(cls, env, default=True):
        """ Downloads and unpacks ADFRsuite. Instructions: -d is the destination folder,
        -c 0/1 picks .pyc/.pyo compilation - no interactive prompt is documented, so no
        'echo "Y" |' is needed. Known issue: on Linux with an older GCC, _openbabel.so can fail with
        "GLIBCXX_3.4.15' not found" - fix is renaming <installDir>/lib/libstdc++.so.6.orig
        to libstdc++.so.6. """
        installer = InstallHelper(ADFRSUITE_DIC['name'], packageHome=cls.getVar(ADFRSUITE_DIC['home']),
                                  packageVersion=ADFRSUITE_DIC['version'])
        installer.addCommand(
            # install.sh needs cwd inside the extracted folder (doesn't cd there itself)
            # running it from outside broke the sibling tarball lookup 
            # (Python2.7.tar.gz etc.) one level up.
            'wget -q https://ccsb.scripps.edu/adfr/download/1038/ -O adfrsuite.tar.gz && '
            'tar -xzf adfrsuite.tar.gz && '
            '(cd ADFRsuite_x86_64Linux_1.0 && ./install.sh -d .. -c 0)',
            targetName=f"{ADFRSUITE_DIC['name']}_installed")
        installer.addPackage(env, dependencies=['wget', 'tar'], default=default)

    @classmethod
    def addVinaPackage(cls, env, default=True):
        """ Installs the Vina CLI binary via conda-forge. PROTAC-Model shells out to
        $VINA/bin/vina, so we need the compiled binary, not just the 'vina' PyPI package
        (Python bindings only). pythonVersion pinned to 3.10, not 3.11: checked
        conda-forge's own repodata. vina=1.2.2 only ships py37-py310 builds, 
        no py311 one, so a 3.11 env would force conda to downgrade the env's own
        Python to satisfy the vina constraint (or fail outright), not the up-front pin
        InstallHelper's own naming implies. """
        installer = InstallHelper(VINA_DIC['name'], packageHome=cls.getVar(VINA_DIC['home']),
                                  packageVersion=VINA_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=VINA_DIC['name'], binaryVersion=VINA_DIC['version'], pythonVersion='3.10'
        ).addCommand(
            f"{cls.getEnvActivationCommand(VINA_DIC)} && conda install -y -c conda-forge vina={VINA_DIC['version']}",
            targetName=f"{VINA_DIC['name']}_installed"
        ).addCommand(
            # getCondaEnvCommand installs under conda's own envs dir, not packageHome
            # without this symlink getVinaProgram() would find an empty folder.
            f"{cls.getEnvActivationCommand(VINA_DIC)} && rm -rf {cls.getVar(VINA_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(VINA_DIC['home'])}",
            targetName=f"{VINA_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addVoromqaPackage(cls, env, default=True):
        """ Installs Voromqa (voronota) via bioconda. """
        installer = InstallHelper(VOROMQA_DIC['name'], packageHome=cls.getVar(VOROMQA_DIC['home']),
                                  packageVersion=VOROMQA_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=VOROMQA_DIC['name'], binaryVersion=VOROMQA_DIC['version'], pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(VOROMQA_DIC)} && conda install -y -c bioconda voronota={VOROMQA_DIC['version']}",
            targetName=f"{VOROMQA_DIC['name']}_installed"
        ).addCommand(
            # Same symlink fix as addVinaPackage above.
            f"{cls.getEnvActivationCommand(VOROMQA_DIC)} && rm -rf {cls.getVar(VOROMQA_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(VOROMQA_DIC['home'])}",
            targetName=f"{VOROMQA_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addFCCPackage(cls, env, default=True):
        """ Clones haddocking/fcc (no PyPI package) and compiles its C/C++ contact
        programs with make. Its Python scripts don't get their own conda env: they're
        only ever invoked from inside PROTAC-Model's own code (pre.cluster()), which
        already runs under PROTAC_MODEL_PYTHON_HOME and passes that interpreter down. """
        installer = InstallHelper(FCC_DIC['name'], packageHome=cls.getVar(FCC_DIC['home']),
                                packageVersion=FCC_DIC['version'])
        installer.getCloneCommand(
            'https://github.com/haddocking/fcc.git', binaryFolderName=FCC_DIC['name'],
            targeName=f"{FCC_DIC['name']}_cloned"
        # FCC's README warns the Makefile may need manual edits, untested.
        ).addCommand(f"cd {FCC_DIC['name']}/src && make", targetName=f"{FCC_DIC['name']}_built")
        installer.addPackage(env, dependencies=['git', 'make', 'gcc'], default=default)

    @classmethod
    def addProtacModelPackage(cls, env, default=True):
        """ Clones gaoqiweng/PROTAC-Model itself, public repo, no build step, its code
        is called directly by the protocol/rosetta/scripts/run_protac_model.py driver. """
        installer = InstallHelper(PROTAC_MODEL_DIC['name'], packageHome=cls.getVar(PROTAC_MODEL_DIC['home']),
                                packageVersion=PROTAC_MODEL_DIC['version'])
        installer.getCloneCommand(
            'https://github.com/gaoqiweng/PROTAC-Model.git', binaryFolderName=PROTAC_MODEL_DIC['name'],
            targeName=f"{PROTAC_MODEL_DIC['name']}_cloned")
        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProtacModelPythonPackage(cls, env, default=True):
        """ Dedicated Python 2.7 conda env (with RDKit) to run PROTAC-Model's own code.
        RDKit's own channel stopped publishing py2.7 builds after 2016.03.3,pinned explicitly below.
        Confirmed on a real install: conda's solver picks a compatible numpy (1.11.3) on
        its own, no manual pin needed. """
        installer = InstallHelper(PROTAC_MODEL_PYTHON_DIC['name'],
                                  packageHome=cls.getVar(PROTAC_MODEL_PYTHON_DIC['home']),
                                  packageVersion=PROTAC_MODEL_PYTHON_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=PROTAC_MODEL_PYTHON_DIC['name'],
            binaryVersion=PROTAC_MODEL_PYTHON_DIC['version'], pythonVersion='2.7'
        ).addCommand(
            f"{cls.getEnvActivationCommand(PROTAC_MODEL_PYTHON_DIC)} && conda install -y -c rdkit rdkit=2016.03.3",
            targetName=f"{PROTAC_MODEL_PYTHON_DIC['name']}_installed"
        ).addCommand(
            # Same symlink fix as addVinaPackage/addVoromqaPackage below.
            f"{cls.getEnvActivationCommand(PROTAC_MODEL_PYTHON_DIC)} && "
            f"rm -rf {cls.getVar(PROTAC_MODEL_PYTHON_DIC['home'])} && "
            f"ln -s $CONDA_PREFIX {cls.getVar(PROTAC_MODEL_PYTHON_DIC['home'])}",
            targetName=f"{PROTAC_MODEL_PYTHON_DIC['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def _requireToolHome(cls, toolDic):
        """ Return toolDic['home']'s configured value, raising if the user never set it.
        Shared by the getXProgram() helpers below so the "not configured" error is
        consistent across tools instead of duplicated once per tool. """
        home = cls.getVar(toolDic['home'])
        if home is None:
            raise FileNotFoundError(
                f"{toolDic['home']} is not set. Point it to your {toolDic['name']} "
                "installation (e.g. in scipion.conf or as a shell environment variable).")
        return home

    # The 4 programs PROTAC-Model's own utils/frodock.py resolves, each shipped as an
    # intel/gcc build pair ('<name>' / '<name>_gcc').
    FRODOCK_BINARIES = ['frodockgrid', 'frodock', 'frodockcluster', 'frodockview']

    @classmethod
    def _binaryLoads(cls, path):
        """ True if path exists and ldd resolves all its shared-library dependencies.
        Existence alone isn't enough: FRODOCK's intel build can be present but missing
        the Intel MKL runtime, which ldd reports as 'not found' rather than failing. """
        if not os.path.exists(path):
            return False
        result = subprocess.run(['ldd', path], capture_output=True, text=True)
        return 'not found' not in result.stdout.lower()

    @classmethod
    def prepareFrodockBinDir(cls, targetDir):
        """ Build a `targetDir/bin/` shim with one symlink per FRODOCK binary, pointing
        each at whichever build (intel or gcc) actually loads here: intel preferred,
        gcc as fallback. PROTAC-Model's own intel/gcc fallback only checks file existence,
        not loadability, so this shim covers that gap. Returns targetDir. """
        home = cls._requireToolHome(FRODOCK_DIC)
        binDir = os.path.join(targetDir, 'bin')
        os.makedirs(binDir, exist_ok=True)

        for name in cls.FRODOCK_BINARIES:
            intel = os.path.join(home, 'bin', name)
            gcc = os.path.join(home, 'bin', f'{name}_gcc')
            if cls._binaryLoads(intel):
                chosen = intel
            # gcc gets the same ldd check, it could be broken too.
            elif cls._binaryLoads(gcc):
                chosen = gcc
            else:
                raise FileNotFoundError(
                    f'Neither {intel} nor {gcc} is usable (checked with ldd).')
            link = os.path.join(binDir, name)
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.abspath(chosen), link)

        # run_protac_model.py copies soap.bin from this shim before running FRODOCK, so
        # it needs to be here too, not just the binaries above.
        soapSrc = os.path.join(home, 'bin', 'soap.bin')
        if not os.path.exists(soapSrc):
            raise FileNotFoundError(f'{soapSrc} not found under FRODOCK_HOME/bin.')
        soapLink = os.path.join(binDir, 'soap.bin')
        if os.path.lexists(soapLink):
            os.remove(soapLink)
        os.symlink(os.path.abspath(soapSrc), soapLink)

        return targetDir

    @classmethod
    def getADFRSuiteProgram(cls, progName):
        """ Return an ADFRsuite binary (reduce, obabel, obenergy, prepare_receptor,
        prepare_ligand...). Unlike FRODOCK, ADFRsuite ships a single build: no intel/gcc
        fallback needed. """
        home = cls._requireToolHome(ADFRSUITE_DIC)
        path = os.path.join(home, 'bin', progName)
        if not os.path.exists(path):
            raise FileNotFoundError(f'{progName} not found under ADFRSUITE_HOME/bin ({home}).')
        return path

    @classmethod
    def getVinaProgram(cls):
        """ Return the Vina binary. No progName parameter: VINA_HOME only ever provides
        this one program, unlike ADFRsuite/FRODOCK which bundle several. """
        home = cls._requireToolHome(VINA_DIC)
        path = os.path.join(home, 'bin', 'vina')
        if not os.path.exists(path):
            raise FileNotFoundError(f'vina not found under VINA_HOME/bin ({home}).')
        return path

    @classmethod
    def getVoromqaProgram(cls):
        """ Return the Voromqa binary (voronota-voromqa). Same single-program case as
        getVinaProgram. """
        home = cls._requireToolHome(VOROMQA_DIC)
        path = os.path.join(home, 'bin', 'voronota-voromqa')
        if not os.path.exists(path):
            raise FileNotFoundError(f'voronota-voromqa not found under VOROMQA_HOME/bin ({home}).')
        return path

    @classmethod
    def getFCCScript(cls, scriptName):
        """ Path to an FCC clustering script (make_contacts.py, calc_fcc_matrix.py,
        cluster_fcc.py, ppretty_clusters.py...). Plain Python 2 scripts, not compiled
        binaries: the caller must still prepend its own interpreter. """
        # FCC_HOME is the parent InstallHelper cloned into. The repo itself lives one
        # level down, in a named subfolder (see addFCCPackage's binaryFolderName).
        home = os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name'])
        path = os.path.join(home, 'scripts', scriptName)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{scriptName} not found under FCC_HOME/{FCC_DIC['name']}/scripts ({home}).")
        return path

    @classmethod
    def getProtacModelScript(cls, scriptName=''):
        """ Path inside the PROTAC-Model checkout, e.g. getProtacModelScript('main.py').
        No scriptName returns the repo root itself. """
        home = os.path.join(cls._requireToolHome(PROTAC_MODEL_DIC), PROTAC_MODEL_DIC['name'])
        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getProtacModelPython(cls):
        """ Path to PROTAC-Model's dedicated Python 2.7+RDKit interpreter. Only used to
        check the env is installed. runCondaScript() activates the env instead of
        launching this path directly. """
        home = cls._requireToolHome(PROTAC_MODEL_PYTHON_DIC)
        path = os.path.join(home, 'bin', 'python2')
        if not os.path.exists(path):
            raise FileNotFoundError(f'python2 not found under PROTAC_MODEL_PYTHON_HOME/bin ({home}).')
        return path

    @classmethod
    def getProtacModelEnviron(cls, frodockHome=None):
        """ Translate our *_HOME variables into the bare names (FRODOCK, VINA, ...)
        PROTAC-Model's own code reads from os.environ, as extraEnvDict for
        runCondaScript(). ROSETTA is required even for a frodock-only run (its module is
        imported unconditionally) and comes from RosettaPlugin, not our own vars.
        frodockHome: pass prepareFrodockBinDir()'s shim to vet intel/gcc; omit for a raw
        FRODOCK_HOME (e.g. from _validate(), which just checks the var is set).
        All paths are made absolute here, since the driver runs from extra/frodock/ or
        extra/rosetta/, not the Scipion project directory a relative path would assume. """
        rosettaHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])
        if rosettaHome is None:
            raise FileNotFoundError(
                f"{ROSETTA_DIC['home']} is not set. Point it to your Rosetta "
                "installation (e.g. in scipion.conf or as a shell environment variable).")
        environ = {
            'FRODOCK': frodockHome or cls._requireToolHome(FRODOCK_DIC),
            'ADFRSUITE': cls._requireToolHome(ADFRSUITE_DIC),
            'VINA': cls._requireToolHome(VINA_DIC),
            'VOROMQA': cls._requireToolHome(VOROMQA_DIC),
            # Same one-level-down layout as getFCCScript().
            'FCC': os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name']),
            'ROSETTA': rosettaHome,
            'PROTAC_MODEL_HOME': cls.getProtacModelScript(),
        }
        return {name: os.path.abspath(path) for name, path in environ.items()}

    @classmethod
    def getPluginScript(cls, scriptName):
        """ Path to a script bundled with this plugin (protac/scripts/<scriptName>).
        Absolute, since it's handed to a process running from another working dir. """
        return os.path.join(_PLUGIN_DIR, 'scripts', scriptName)
    
    @classmethod
    def runCondaScript(cls, scriptPath, args, condaDic, extraEnvDict=None, cwd=None):
        """ Run a Python script with condaDic's env activated first, rather than calling
        that env's interpreter by absolute path - needed because PROTAC-Model itself
        shells out to FCC's scripts as a bare 'python ...', which only finds the right
        Python 2.7 if the env's bin/ is on PATH. scriptPath is made absolute here since
        cwd is not the Scipion project directory. """
        program = f'{cls.getEnvActivationCommand(condaDic)} && python "{os.path.abspath(scriptPath)}"'
        cls.runProgram(program, args, extraEnvDict=extraEnvDict, cwd=cwd)

    @classmethod
    def getEnviron(cls):
        """ Base environment for launching external programs. Tool-specific vars are
        layered on top via runProgram's extraEnvDict, not here. """
        return pwutils.Environ(os.environ)

    @classmethod
    def runProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Launch an external program with the given env/cwd. Not tool-specific.
        cwd is resolved to absolute here so failures log an unambiguous path. """
        env = cls.getEnviron()
        if extraEnvDict is not None:
            env.update(extraEnvDict)
        pwutils.runJob(None, program, args, env=env,
                       cwd=os.path.abspath(cwd) if cwd else cwd)
