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

import pyworkflow.utils as pwutils
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from rosetta import Plugin as RosettaPlugin, ROSETTA_DIC

from .constants import *

_version_ = "0.1"
# FRODOCK is a separate external tool (not Rosetta), used by the PROTAC-Model pipeline
# for the initial global protein-protein docking step. No 'version' key: we don't pin/
# validate a specific FRODOCK version, only that FRODOCK_HOME points somewhere real.
FRODOCK_DIC = {'name': 'frodock', 'home': 'FRODOCK_HOME'}

# The remaining four are used later, by filterPosesStep (PROTAC-Model's filter_frodock()):
# reduce/obabel/obenergy/prepare_receptor/prepare_ligand (ADFRsuite), vina (Vina), the
# voronota-voromqa binary (Voromqa), and the FCC clustering scripts. Unlike FRODOCK/Rosetta
# (license-gated downloads, see below), none of these four require a license click-through,
# so they get real defineBinaries() support via InstallHelper - 'version' is needed now
# (InstallHelper/getEnvName use it to name the conda env / package folder).
ADFRSUITE_DIC = {'name': 'adfrsuite', 'version': '1.0', 'home': 'ADFRSUITE_HOME'}
VINA_DIC = {'name': 'vina', 'version': '1.2.7', 'home': 'VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'version': '1.29.4816', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'version': 'latest', 'home': 'FCC_HOME'}

# PROTAC-Model's own code (main.py, utils/*) is called directly, not reimplemented (see
# protocol_protac_model.py/protac/scripts/run_*.py). Its code is genuine Python 2, so
# PROTAC_MODEL_PYTHON_HOME is a dedicated Python 2.7+rdkit conda env, never scipion3's own.
PROTAC_MODEL_DIC = {'name': 'protac-model', 'version': 'latest', 'home': 'PROTAC_MODEL_HOME'}
PROTAC_MODEL_PYTHON_DIC = {'name': 'protac-model-python', 'version': '2.7',
                           'home': 'PROTAC_MODEL_PYTHON_HOME'}

class Plugin(pwchemPlugin):
    _homeVar = PROTAC_MODEL_DIC['home']
    _pathVars = [PROTAC_MODEL_DIC['home']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. """
        # FRODOCK_HOME stays manual (license click-through, see FRODOCK_DIC above) - None
        # until the user points it at their own download.
        cls._defineVar(FRODOCK_DIC['home'], None)

        # The five tools below have no license gate - _defineEmVar wires each home to
        # wherever defineBinaries()/InstallHelper installs it.
        cls._defineEmVar(ADFRSUITE_DIC['home'], cls.getEnvName(ADFRSUITE_DIC))
        cls._defineEmVar(VINA_DIC['home'], cls.getEnvName(VINA_DIC))
        cls._defineEmVar(VOROMQA_DIC['home'], cls.getEnvName(VOROMQA_DIC))
        cls._defineEmVar(FCC_DIC['home'], cls.getEnvName(FCC_DIC))
        cls._defineEmVar(PROTAC_MODEL_DIC['home'], cls.getEnvName(PROTAC_MODEL_DIC))
        cls._defineEmVar(PROTAC_MODEL_PYTHON_DIC['home'], cls.getEnvName(PROTAC_MODEL_PYTHON_DIC))

    @classmethod
    def defineBinaries(cls, env):
        # FRODOCK/Rosetta excluded on purpose (manual, license-gated - see FRODOCK_DIC
        # above; ROSETTA_HOME is scipion-chem-rosetta's own responsibility).
        cls.addADFRSuitePackage(env)
        cls.addVinaPackage(env)
        cls.addVoromqaPackage(env)
        cls.addFCCPackage(env)
        cls.addProtacModelPackage(env)
        cls.addProtacModelPythonPackage(env)

    # ---------------------------- Package installers (InstallHelper) -------------
    @classmethod
    def addADFRSuitePackage(cls, env, default=True):
        """ Downloads and unpacks ADFRsuite (no license click-through found on
        ccsb.scripps.edu for the non-commercial installer, unlike FRODOCK).
        URL and install.sh flags confirmed against ccsb.scripps.edu/adfr/downloads/'s own
        "INSTALLING FROM TARBALL" instructions (2026-09-10): -d is the destination folder,
        -c 0/1 picks .pyc/.pyo compilation - no interactive prompt is documented, so no
        'echo "Y" |' is needed.
        Known issue (same page): on Linux with an older GCC, _openbabel.so can fail with
        "GLIBCXX_3.4.15' not found" - fix is renaming <installDir>/lib/libstdc++.so.6.orig
        to libstdc++.so.6. """
        installer = InstallHelper(ADFRSUITE_DIC['name'], packageHome=cls.getVar(ADFRSUITE_DIC['home']),
                                  packageVersion=ADFRSUITE_DIC['version'])
        installer.addCommand(
            # install.sh needs cwd inside the extracted folder (doesn't cd there itself) -
            # confirmed on the CNB VM (2026-09-10), running it from outside broke the
            # sibling tarball lookup (Python2.7.tar.gz etc.) one level up.
            'wget -q https://ccsb.scripps.edu/adfr/download/1038/ -O adfrsuite.tar.gz && '
            'tar -xzf adfrsuite.tar.gz && '
            '(cd ADFRsuite_x86_64Linux_1.0 && ./install.sh -d .. -c 0)',
            targetName=f"{ADFRSUITE_DIC['name']}_installed")
        installer.addPackage(env, dependencies=['wget', 'tar'], default=default)

    @classmethod
    def addVinaPackage(cls, env, default=True):
        """ Installs the Vina CLI binary via conda-forge - PROTAC-Model shells out to
        $VINA/bin/vina, so we need the compiled binary, not just the 'vina' PyPI package
        (Python bindings only). """
        installer = InstallHelper(VINA_DIC['name'], packageHome=cls.getVar(VINA_DIC['home']),
                                  packageVersion=VINA_DIC['version'])
        installer.getCondaEnvCommand(
            binaryName=VINA_DIC['name'], binaryVersion=VINA_DIC['version'], pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(VINA_DIC)} && conda install -y -c conda-forge vina={VINA_DIC['version']}",
            targetName=f"{VINA_DIC['name']}_installed"
        ).addCommand(
            # getCondaEnvCommand installs under conda's own envs dir, not packageHome -
            # without this symlink getVinaProgram() would find an empty folder. Same fix
            # pwchem uses for MGLTools (pwchem/__init__.py, addMGLToolsPackage).
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
        # binaryFolderName=FCC_DIC['name'] (not '.'): getCloneCommand's 'cd <packageHome>
        # && git clone <url> .' only works if packageHome already exists and is empty -
        # pwchem always clones into a named subfolder instead (see e.g. addShapeItPackage).
        installer.getCloneCommand(
            'https://github.com/haddocking/fcc.git', binaryFolderName=FCC_DIC['name'],
            targeName=f"{FCC_DIC['name']}_cloned"
        # XXX FCC's README warns the Makefile may need manual edits - untested.
        ).addCommand(f"cd {FCC_DIC['name']}/src && make", targetName=f"{FCC_DIC['name']}_built")
        installer.addPackage(env, dependencies=['git', 'make', 'gcc'], default=default)

    @classmethod
    def addProtacModelPackage(cls, env, default=True):
        """ Clones gaoqiweng/PROTAC-Model itself - public repo, no build step, its code
        is called directly by the protocol/rosetta/scripts/run_*.py drivers. """
        installer = InstallHelper(PROTAC_MODEL_DIC['name'], packageHome=cls.getVar(PROTAC_MODEL_DIC['home']),
                                  packageVersion=PROTAC_MODEL_DIC['version'])
        # Same binaryFolderName reasoning as addFCCPackage above (named subfolder, not
        # '.') - getProtacModelScript()/getProtacModelPython() below account for it.
        installer.getCloneCommand(
            'https://github.com/gaoqiweng/PROTAC-Model.git', binaryFolderName=PROTAC_MODEL_DIC['name'],
            targeName=f"{PROTAC_MODEL_DIC['name']}_cloned")
        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProtacModelPythonPackage(cls, env, default=True):
        """ Dedicated Python 2.7 conda env (with RDKit) to run PROTAC-Model's own code -
        never scipion3's own env, and not the local 'protac-model' conda env some
        machines may already have lying around (that one is Python 3.10, incompatible
        with this genuinely-Python-2 codebase).
        RDKit's own channel stopped publishing py2.7 builds after 2016.03.3 (checked
        anaconda.org/rdkit/rdkit's full file list, 2026-09-10) - pinned explicitly below.
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

    @classmethod
    def getFrodockProgram(cls, progName):
        """ Return the FRODOCK binary that will be used, trying the intel build first and
        falling back to the gcc build (FRODOCK ships both, e.g. frodockgrid/frodockgrid_gcc,
        and only one is guaranteed to work on a given machine). """
        home = cls._requireToolHome(FRODOCK_DIC)

        # The two candidate paths, same convention as FRODOCK's own <name>/<name>_gcc pair.
        intel = os.path.join(home, 'bin', progName)
        gcc = os.path.join(home, 'bin', f'{progName}_gcc')
        if os.path.exists(intel):
            return intel
        elif os.path.exists(gcc):
            return gcc
        else:
            # Fail loudly with both checked paths, instead of the original PROTAC-Model
            # script's print() + sys.exit() (which would kill the whole Scipion process).
            raise FileNotFoundError(
                f'{progName} not found under FRODOCK_HOME/bin ({home}). Checked {intel} and {gcc}.')

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
        # FCC_HOME is the parent InstallHelper cloned into; the repo itself lives one
        # level down, in a named subfolder (see addFCCPackage's binaryFolderName).
        home = os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name'])
        path = os.path.join(home, 'scripts', scriptName)
        if not os.path.exists(path):
            raise FileNotFoundError(f"{scriptName} not found under FCC_HOME/{FCC_DIC['name']}/scripts ({home}).")
        return path

    @classmethod
    def getProtacModelScript(cls, scriptName=''):
        """ Return a path inside the PROTAC-Model checkout (gaoqiweng/PROTAC-Model), e.g.
        getProtacModelScript() for the repo root (what our own rosetta/scripts/
        run_protac_model.py driver adds to sys.path), or getProtacModelScript('main.py')
        for a specific file. Same subfolder reasoning as getFCCScript above: PROTAC_MODEL_HOME
        is the parent InstallHelper cloned into, the repo itself is one level down. """
        home = os.path.join(cls._requireToolHome(PROTAC_MODEL_DIC), PROTAC_MODEL_DIC['name'])
        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getProtacModelPython(cls):
        """ Path to the Python 2.7 (+ RDKit) interpreter dedicated to PROTAC-Model's own
        code - only used to check the env is installed (see _validate in the protocol),
        not to launch anything: runCondaScript() below activates this same conda env
        instead, so the env's own 'python' lands on PATH (a bare absolute-interpreter
        launch wouldn't, and FCC's clustering scripts need that). """
        home = cls._requireToolHome(PROTAC_MODEL_PYTHON_DIC)
        path = os.path.join(home, 'bin', 'python2')
        if not os.path.exists(path):
            raise FileNotFoundError(f'python2 not found under PROTAC_MODEL_PYTHON_HOME/bin ({home}).')
        return path

    @classmethod
    def getProtacModelEnviron(cls):
        """ PROTAC-Model's own utils/*.py modules read these bare names (no _HOME suffix)
        from os.environ at import time - this translates our *_HOME variables into that
        convention in one place, for use as runCondaScript()'s extraEnvDict. ROSETTA is
        required unconditionally (not just when refining): utils/rosetta.py is imported
        by our driver script at module load for every phase, so ROSETTA must resolve even
        for a frodock-only run. Unlike the other five, ROSETTA_HOME is owned by the
        scipion-chem-rosetta plugin (a real dependency, see requirements.txt), so it's
        resolved via RosettaPlugin.getVar() instead of our own _requireToolHome(). """
        rosettaHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])
        if rosettaHome is None:
            raise FileNotFoundError(
                f"{ROSETTA_DIC['home']} is not set. Point it to your Rosetta "
                "installation (e.g. in scipion.conf or as a shell environment variable).")
        return {
            'FRODOCK': cls._requireToolHome(FRODOCK_DIC),
            'ADFRSUITE': cls._requireToolHome(ADFRSUITE_DIC),
            'VINA': cls._requireToolHome(VINA_DIC),
            'VOROMQA': cls._requireToolHome(VOROMQA_DIC),
            # Same one-level-down layout as getFCCScript().
            'FCC': os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name']),
            'ROSETTA': rosettaHome,
            'PROTAC_MODEL_HOME': cls.getProtacModelScript(),
        }

    @classmethod
    def getPluginScript(cls, scriptName):
        """ Path to a script bundled with this plugin itself (protac/scripts/<scriptName>),
        e.g. run_protac_model.py - mirrors pwchem's Plugin.getScriptsDir(). """
        return os.path.join(os.path.dirname(__file__), 'scripts', scriptName)
    @classmethod
    def runCondaScript(cls, scriptPath, args, condaDic, extraEnvDict=None, cwd=None):
        """ Launch a Python script with a conda env activated first, instead of calling
        that env's interpreter by absolute path - mirrors pwchem's Plugin.runScript(). The
        difference matters here: run_protac_model.py (--phase filter) ends up shelling out
        to FCC's clustering scripts as a bare 'python <script>.py' command (see
        PROTAC-Model's own preprocess.py) - that only resolves to the right Python 2.7 if
        this env's bin/ is actually on PATH, which activating it does and launching by
        absolute interpreter path alone would not. """
        program = f'{cls.getEnvActivationCommand(condaDic)} && python {scriptPath}'
        cls.runProgram(program, args, extraEnvDict=extraEnvDict, cwd=cwd)

    @classmethod
    def getEnviron(cls):
        """ Base environment for launching external programs - starts from the current
        process environment; tool-specific additions (FRODOCK/ADFRSUITE/VINA/...) are
        layered on top via runProgram's extraEnvDict, not here. """
        return pwutils.Environ(os.environ)

    @classmethod
    def runProgram(cls, program, args=None, extraEnvDict=None, cwd=None):
        """ Internal shortcut function to launch an external program (Rosetta or, e.g.,
        FRODOCK). Not tool-specific: only builds the environment and launches the process. """
        env = cls.getEnviron()
        if extraEnvDict is not None:
            env.update(extraEnvDict)
        pwutils.runJob(None, program, args, env=env, cwd=cwd)
