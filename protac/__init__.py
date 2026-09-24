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

import glob
import os
import subprocess

import pyworkflow.utils as pwutils
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from pwchem import constants as pwchemConstants
from rosetta import Plugin as RosettaPlugin, ROSETTA_DIC

from .constants import *

# Computed once at import, before any chdir() the protocol's steps do later,
# abspath(__file__) could resolve wrong if computed after the cwd has moved.
_PLUGIN_DIR = os.path.dirname(os.path.abspath(__file__))

__version__ = ALPHA_VERSION


class Plugin(pwchemPlugin):
    _homeVar = PROTAC_MODEL_DIC['home']
    _pathVars = [PROTAC_MODEL_DIC['home']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file. """
        # Defaults point where defineBinaries() installs each tool; a value set in
        # scipion.conf or the environment still wins.
        cls._defineEmVar(FRODOCK_DIC['home'], cls.getEnvName(FRODOCK_DIC))
        cls._defineEmVar(ADFRSUITE_DIC['home'], cls.getEnvName(ADFRSUITE_DIC))
        cls._defineEmVar(VINA_DIC['home'], cls.getEnvName(VINA_DIC))
        cls._defineEmVar(VOROMQA_DIC['home'], cls.getEnvName(VOROMQA_DIC))
        cls._defineEmVar(FCC_DIC['home'], cls.getEnvName(FCC_DIC))
        cls._defineEmVar(PROTAC_MODEL_DIC['home'], cls.getEnvName(PROTAC_MODEL_DIC))
        cls._defineEmVar(PROTAC_MODEL_PYTHON_DIC['home'], cls.getEnvName(PROTAC_MODEL_PYTHON_DIC))
        cls._defineEmVar(PROSETTAC_DIC['home'], cls.getEnvName(PROSETTAC_DIC))
        cls._defineEmVar(PROSETTAC_PYTHON_DIC['home'], cls.getEnvName(PROSETTAC_PYTHON_DIC))
        cls._defineEmVar(PROSETTAC_PYTHON2_DIC['home'], cls.getEnvName(PROSETTAC_PYTHON2_DIC))
        # _defineVar, not _defineEmVar: must be able to stay None so
        # _requirePatchdockHome() can tell "not installed" apart from "installed".
        cls._defineVar(PATCHDOCK_DIC['home'], None)

    @classmethod
    def defineBinaries(cls, env):
        # Rosetta and PatchDock excluded on purpose: both need a personal academic
        # registration. ROSETTA_HOME is scipion-chem-rosetta's own responsibility.
        cls.addFrodockPackage(env)
        cls.addADFRSuitePackage(env)
        cls.addVinaPackage(env)
        cls.addVoromqaPackage(env)
        cls.addFCCPackage(env)
        cls.addProtacModelPackage(env)
        cls.addProtacModelPythonPackage(env)
        cls.addProsettaCPackage(env)
        cls.addProsettaCPythonPackage(env)
        cls.addProsettaCPython2Package(env)

    # ---------------------------- Package installers (InstallHelper) -------------
    @classmethod
    def addFrodockPackage(cls, env, default=True):
        """ Downloads FRODOCK's prebuilt intel and gcc binaries. The URL carries a token
        that changes when the vendor updates the site, and a stale one returns HTTP 200
        with an HTML page, hence the gzip check. """
        installer = InstallHelper(FRODOCK_DIC['name'], packageHome=cls.getVar(FRODOCK_DIC['home']),
                                  packageVersion=FRODOCK_DIC['version'])
        installer.addCommand(
            # Single-quoted: the URL carries '&' and '[]', and these commands are run
            # through a shell, so unquoted the '&' would split it into background jobs.
            "wget -q -O frodock.tgz 'https://chaconlab.org/component/zoo/"
            "?task=callelement&format=raw&item_id=17"
            "&element=f85c494b-2b32-4109-b8c1-083cca2b7db6&method=download"
            "&args[0]=3acfd359e585affd517dbfe436236163' && "
            # The exit below kills this whole shell outright on a bad download, so no
            # outer grouping is needed to keep 'tar'/'rm' from running afterwards - only
            # the inner '{ }' is needed, to make echo+exit a single branch of the '||'.
            "gzip -t frodock.tgz 2>/dev/null || { echo 'FRODOCK download failed: the "
            "server did not return a tarball. The download link carries a token that "
            "changes when the vendor updates their site; get the current one by copying "
            "the target of the Download button at "
            "https://chaconlab.org/modeling/frodock/frodock-donwload and update it in "
            "addFrodockPackage, or install FRODOCK by hand and point FRODOCK_HOME at "
            "it.' >&2; exit 1; } && "
            # --strip-components=1 drops the tarball's own frodock3_linux64/ top level, so
            # the binaries land in <FRODOCK_HOME>/bin, which is where the rest of this
            # plugin (prepareFrodockBinDir) looks for them.
            "tar -xzf frodock.tgz --strip-components=1 && "
            "rm -f frodock.tgz",
            targetName=f"{FRODOCK_DIC['name']}_installed")
        installer.addPackage(env, dependencies=['wget', 'tar', 'gzip'], default=default)

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
    def _addCondaEnvPackage(cls, env, dic, pythonVersion, installCmd=None, default=True):
        """ Creates a conda env for dic, optionally runs installCmd inside it, and symlinks
        the env into dic's home: conda puts envs under its own envs/ dir, while the rest
        of the plugin looks for them at the home _defineEmVar composes. """
        activation = cls.getEnvActivationCommand(dic)
        home = cls.getVar(dic['home'])
        installer = InstallHelper(dic['name'], packageHome=home, packageVersion=dic['version'])
        installer.getCondaEnvCommand(binaryName=dic['name'], binaryVersion=dic['version'],
                                     pythonVersion=pythonVersion)
        if installCmd:
            installer.addCommand(f"{activation} && {installCmd}",
                                 targetName=f"{dic['name']}_installed")
        installer.addCommand(f"{activation} && rm -rf {home} && ln -s $CONDA_PREFIX {home}",
                             targetName=f"{dic['name']}_symlinked")
        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addVinaPackage(cls, env, default=True):
        """ The Vina CLI binary from conda-forge: PROTAC-Model shells out to
        $VINA/bin/vina, which the 'vina' PyPI package (bindings only) doesn't provide. """
        cls._addCondaEnvPackage(env, VINA_DIC, '3.10',
                                f"conda install -y -c conda-forge vina={VINA_DIC['version']}",
                                default=default)

    @classmethod
    def addVoromqaPackage(cls, env, default=True):
        cls._addCondaEnvPackage(env, VOROMQA_DIC, '3.11',
                                f"conda install -y -c bioconda voronota={VOROMQA_DIC['version']}",
                                default=default)

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
        """ Python 2.7 env with RDKit for PROTAC-Model's own code. 2016.03.3 is the last
        py2.7 build on RDKit's channel; conda picks a compatible numpy by itself. """
        cls._addCondaEnvPackage(env, PROTAC_MODEL_PYTHON_DIC, '2.7',
                                'conda install -y -c rdkit rdkit=2016.03.3', default=default)

    @classmethod
    def addProsettaCPackage(cls, env, default=True):
        """ Clones LondonLab/PRosettaC itself, public repo, no build step. """
        installer = InstallHelper(PROSETTAC_DIC['name'], packageHome=cls.getVar(PROSETTAC_DIC['home']),
                                  packageVersion=PROSETTAC_DIC['version'])
        installer.getCloneCommand(
            'https://github.com/LondonLab/PRosettaC.git', binaryFolderName=PROSETTAC_DIC['name'],
            targeName=f"{PROSETTAC_DIC['name']}_cloned")
        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProsettaCPythonPackage(cls, env, default=True):
        """ Py3 env for PRosettaC's own code. RDKit 2023.09 is built against numpy 1.x. """
        cls._addCondaEnvPackage(env, PROSETTAC_PYTHON_DIC, '3.10',
                                'conda install -y -c conda-forge rdkit=2023.09.6 '
                                'numpy=1.26.4 scikit-learn=1.5.2',
                                default=default)

    @classmethod
    def addProsettaCPython2Package(cls, env, default=True):
        """ Bare Python 2.7 env: PRosettaC runs molfile_to_params.py as a literal
        'python2.7 ...', so this env's bin/ goes on PATH (see getPRosettaCEnviron). """
        cls._addCondaEnvPackage(env, PROSETTAC_PYTHON2_DIC, '2.7', default=default)

    @classmethod
    def _requireToolHome(cls, toolDic, binary=''):
        """ toolDic's home, raising if it (or binary inside it) is missing on disk. """
        home = cls.getVar(toolDic['home'])
        if home is None or not os.path.exists(os.path.join(home, binary)):
            raise FileNotFoundError(
                f"{toolDic['home']} does not point to an existing {toolDic['name']} "
                f"installation (got: {home}). Install it with 'scipion3 installb "
                f"{cls.getEnvName(toolDic)}', or set the variable (e.g. in scipion.conf "
                "or as a shell environment variable) to your own installation.")
        return home

    @classmethod
    def _requirePatchdockHome(cls):
        """ Like _requireToolHome, but PatchDock isn't scipion3-installable, so the error
        points at manual install instead of 'scipion3 installb'. """
        home = cls.getVar(PATCHDOCK_DIC['home'])
        if home is None or not os.path.isdir(home):
            raise FileNotFoundError(
                f"{PATCHDOCK_DIC['home']} is not set. PatchDock needs a separate "
                "academic license/registration (https://bioinfo3d.cs.tau.ac.il/PatchDock/) "
                f"and cannot be installed automatically. Install it yourself and point "
                f"{PATCHDOCK_DIC['home']} at it. Got: {home}.")
        return home

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

        for name in FRODOCK_BINARIES:
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
            cls._forceSymlink(chosen, os.path.join(binDir, name))

        # run_protac_model.py copies soap.bin from this shim before running FRODOCK, so
        # it needs to be here too, not just the binaries above.
        soapSrc = os.path.join(home, 'bin', 'soap.bin')
        if not os.path.exists(soapSrc):
            raise FileNotFoundError(f'{soapSrc} not found under FRODOCK_HOME/bin.')
        cls._forceSymlink(soapSrc, os.path.join(binDir, 'soap.bin'))

        return targetDir

    @classmethod
    def _forceSymlink(cls, source, link):
        """ Symlink source at link, replacing whatever was there (file, link or a real
        directory left by an older shim layout). """
        pwutils.makePath(os.path.dirname(link))
        pwutils.cleanPath(link)
        pwutils.createAbsLink(source, link)
        return link

    @classmethod
    def requireRosettaScriptsBinary(cls):
        """ A rosetta_scripts binary under ROSETTA_HOME. PRosettaC hardcodes the .default.
        build, but installs often ship only .static. or .mpi., so any build is accepted
        and the shim renames it. """
        home = RosettaPlugin.getVar(ROSETTA_DIC['home'])
        if home is None or not os.path.isdir(home):
            raise FileNotFoundError(
                f"{ROSETTA_DIC['home']} does not point to an existing Rosetta "
                f'installation (got: {home}). Set it (e.g. in scipion.conf or as a shell '
                'environment variable) to your own installation.')

        binDir = os.path.join(home, 'main', 'source', 'bin')
        for build in ROSETTA_SCRIPTS_BUILDS:
            path = os.path.join(binDir, f'rosetta_scripts.{build}.linuxgccrelease')
            if os.path.exists(path):
                return path
        # Any other build will do; sorted so the choice is reproducible.
        others = sorted(glob.glob(os.path.join(binDir, 'rosetta_scripts.*')))
        if others:
            return others[0]
        raise FileNotFoundError(
            f'No rosetta_scripts binary found in {binDir}. PRosettaC needs a compiled '
            'Rosetta, not just the source bundle.')

    @classmethod
    def prepareRosettaScriptsShim(cls, targetDir):
        """ Builds targetDir as a stand-in ROSETTA_FOL: symlinks into the real install
        with the layout PRosettaC expects, exposing the binary under the .default. name.
        The scripts, tools and database are linked too, since PRosettaC also resolves
        molfile_to_params.py and clean_pdb.py against ROSETTA_FOL. """
        realHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])
        cls._forceSymlink(cls.requireRosettaScriptsBinary(), os.path.join(
            targetDir, 'main', 'source', 'bin', 'rosetta_scripts.default.linuxgccrelease'))

        # Whole directories, so scripts that import their siblings still work.
        for relPath in (os.path.join('main', 'source', 'scripts'),
                        os.path.join('main', 'database'), 'tools'):
            source = os.path.join(realHome, relPath)
            if not os.path.isdir(source):
                raise FileNotFoundError(
                    f'{source} not found under ROSETTA_HOME. PRosettaC needs the full '
                    'Rosetta bundle (binaries, python scripts, tools and database), not '
                    'just the binaries.')
            cls._forceSymlink(source, os.path.join(targetDir, relPath))
        return targetDir

    @classmethod
    def requireObabelBinary(cls):
        """ pwchem installs OpenBabel as a named conda env, so OPENBABEL_HOME has no bin/:
        resolve the binary inside the env instead. """
        obabelBin = pwchemPlugin.getEnvPath(pwchemConstants.OPENBABEL_DIC, innerPath='bin/obabel')
        if not os.path.exists(obabelBin):
            raise FileNotFoundError(f"obabel not found at {obabelBin}. Is pwchem's "
                                    'OpenBabel environment installed?')
        return obabelBin

    @classmethod
    def preparePRosettaCBabelShim(cls, targetDir):
        """ Builds targetDir/babel, translating PRosettaC's OpenBabel 2.x calls
        ('babel in out [-h]') into pwchem's obabel 3 ('obabel in -O out [-h]'). """
        obabelBin = cls.requireObabelBinary()
        os.makedirs(targetDir, exist_ok=True)
        wrapperPath = os.path.join(targetDir, 'babel')
        activate = pwchemPlugin.getEnvActivationCommand(pwchemConstants.OPENBABEL_DIC)
        script = (
            '#!/usr/bin/env bash\n'
            # No '-u'/pipefail: conda's activation scripts would trip them.
            'set -e\n'
            # obabel finds its plugins through the env; the absolute path avoids picking
            # another obabel from PATH.
            f'{activate}\n'
            f'OBABEL="{os.path.abspath(obabelBin)}"\n'
            'IN="$1"; OUT="$2"; shift 2 || true\n'
            # in==out (round-trip call) needs a temp file, can't overwrite while reading.
            'if [ "$IN" = "$OUT" ]; then\n'
            '  TMP="$(mktemp --suffix=".${OUT##*.}")"\n'
            '  "$OBABEL" "$IN" -O "$TMP" "$@"\n'
            '  mv "$TMP" "$OUT"\n'
            'else\n'
            '  "$OBABEL" "$IN" -O "$OUT" "$@"\n'
            'fi\n')
        with open(wrapperPath, 'w') as f:
            f.write(script)
        os.chmod(wrapperPath, 0o755)
        return targetDir

    @classmethod
    def getProtacModelScript(cls, scriptName=''):
        """ Path inside the PROTAC-Model checkout, e.g. getProtacModelScript('main.py').
        No scriptName returns the repo root itself. """
        home = os.path.join(cls._requireToolHome(PROTAC_MODEL_DIC), PROTAC_MODEL_DIC['name'])
        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getProtacModelPython(cls):
        """ PROTAC-Model's Python 2.7 interpreter, only used to check the env exists. """
        home = cls._requireToolHome(PROTAC_MODEL_PYTHON_DIC)
        path = os.path.join(home, 'bin', 'python2')
        if not os.path.exists(path):
            raise FileNotFoundError(f'python2 not found under PROTAC_MODEL_PYTHON_HOME/bin ({home}).')
        return path

    @classmethod
    def getProtacModelEnviron(cls, frodockHome=None):
        """ Our *_HOME variables under the names PROTAC-Model reads (FRODOCK, VINA...),
        as absolute paths. ROSETTA is needed even without refinement, because its module
        is always imported. frodockHome: the shim from prepareFrodockBinDir(), if built. """
        rosettaHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])
        if rosettaHome is None:
            raise FileNotFoundError(
                f"{ROSETTA_DIC['home']} is not set. Point it to your Rosetta "
                "installation (e.g. in scipion.conf or as a shell environment variable).")
        environ = {
            'FRODOCK': frodockHome or cls._requireToolHome(FRODOCK_DIC),
            'ADFRSUITE': cls._requireToolHome(ADFRSUITE_DIC),
            'VINA': cls._requireToolHome(VINA_DIC, 'bin/vina'),
            'VOROMQA': cls._requireToolHome(VOROMQA_DIC),
            # The clone lives one level down, in a subfolder named after the package.
            'FCC': os.path.join(cls._requireToolHome(FCC_DIC), FCC_DIC['name']),
            'ROSETTA': rosettaHome,
            'PROTAC_MODEL_HOME': cls.getProtacModelScript(),
        }
        return {name: os.path.abspath(path) for name, path in environ.items()}

    @classmethod
    def getPRosettaCScript(cls, scriptName=''):
        """ Path inside the PRosettaC checkout. No scriptName returns the repo root. """
        home = os.path.join(cls._requireToolHome(PROSETTAC_DIC), PROSETTAC_DIC['name'])
        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getPRosettaCEnviron(cls, rosettaHome, obDir):
        """ The 4 variables PRosettaC reads at import time, plus the Python 2.7 env on PATH
        for its bare 'python2.7' call. SCRIPTS_FOL needs a trailing slash, as rosetta.py
        appends to it directly. rosettaHome/obDir: the two shims. """
        python2Home = cls._requireToolHome(PROSETTAC_PYTHON2_DIC)
        environ = {
            'PATCHDOCK': cls._requirePatchdockHome(),
            'OB': os.path.abspath(obDir),
            'SCRIPTS_FOL': os.path.abspath(cls.getPRosettaCScript()) + os.sep,
            'ROSETTA_FOL': os.path.abspath(rosettaHome),
            'PATH': os.path.join(python2Home, 'bin') + os.pathsep + os.environ.get('PATH', ''),
        }
        return environ

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
