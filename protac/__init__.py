import glob
import os
import subprocess

import pyworkflow.utils as pwutils
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from pwchem import constants as pwchemConstants
from rosetta import Plugin as RosettaPlugin, ROSETTA_DIC

from .constants import *


_PLUGIN_DIR = os.path.dirname(os.path.abspath(__file__))

__version__ = ALPHA_VERSION


class Plugin(pwchemPlugin):

    # Default installation paths. ADFRSuite and Vina are provided by
    # scipion-chem-autodock; this plugin does not install them.
    _frodockHome = pwchemPlugin.getDefPath(FRODOCK_DIC)
    _adfrsuiteHome = pwchemPlugin.getDefPath(ADFRSUITE_DIC)
    _vinaHome = pwchemPlugin.getDefPath(VINA_DIC)
    _voromqaHome = pwchemPlugin.getDefPath(VOROMQA_DIC)
    _fccHome = pwchemPlugin.getDefPath(FCC_DIC)
    _protacModelHome = pwchemPlugin.getDefPath(PROTAC_MODEL_DIC)
    _protacModelPythonHome = pwchemPlugin.getDefPath(PROTAC_MODEL_PYTHON_DIC)
    _prosettaCHome = pwchemPlugin.getDefPath(PROSETTAC_DIC)
    _prosettaCPythonHome = pwchemPlugin.getDefPath(PROSETTAC_PYTHON_DIC)
    _prosettaCPython2Home = pwchemPlugin.getDefPath(PROSETTAC_PYTHON2_DIC)

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(FRODOCK_DIC['home'], cls._frodockHome)
        cls._defineEmVar(ADFRSUITE_DIC['home'], cls._adfrsuiteHome)
        cls._defineEmVar(VINA_DIC['home'], cls._vinaHome)
        cls._defineEmVar(VOROMQA_DIC['home'], cls._voromqaHome)
        cls._defineEmVar(FCC_DIC['home'], cls._fccHome)
        cls._defineEmVar(PROTAC_MODEL_DIC['home'], cls._protacModelHome)
        cls._defineEmVar(
            PROTAC_MODEL_PYTHON_DIC['home'], cls._protacModelPythonHome)
        cls._defineEmVar(PROSETTAC_DIC['home'], cls._prosettaCHome)
        cls._defineEmVar(
            PROSETTAC_PYTHON_DIC['home'], cls._prosettaCPythonHome)
        cls._defineEmVar(
            PROSETTAC_PYTHON2_DIC['home'], cls._prosettaCPython2Home)

        # PatchDock is optional and must be installed manually.
        cls._defineVar(PATCHDOCK_DIC['home'], None)

    @classmethod
    def defineBinaries(cls, env):
        # FRODOCK and PROTAC-specific dependencies are installed here.
        # ADFRSuite and Vina are installed by scipion-chem-autodock.
        # Rosetta and PatchDock require separate installations.
        #cls.addFrodockPackage(env)
        cls.addVoromqaPackage(env)
        cls.addFCCPackage(env)
        cls.addProtacModelPackage(env)
        cls.addProtacModelPythonPackage(env)
        cls.addProsettaCPackage(env)
        cls.addProsettaCPythonPackage(env)
        cls.addProsettaCPython2Package(env)

    # -------------------------------------------------------------------------
    # Package installers
    # -------------------------------------------------------------------------

    @classmethod
    def addFrodockPackage(cls, env, default=True):
        installer = InstallHelper(
            FRODOCK_DIC['name'],
            packageHome=cls.getVar(FRODOCK_DIC['home']),
            packageVersion=FRODOCK_DIC['version'])

        installer.addCommand(
            "wget -q -O frodock.tgz "
            "'https://chaconlab.org/component/zoo/"
            "?task=callelement&format=raw&item_id=17"
            "&element=f85c494b-2b32-4109-b8c1-083cca2b7db6"
            "&method=download&args[0]=3acfd359e585affd517dbfe436236163' "
            "&& gzip -t frodock.tgz 2>/dev/null || "
            "{ echo 'FRODOCK download failed: check the current vendor "
            "download URL or install FRODOCK manually and set FRODOCK_HOME.' "
            ">&2; exit 1; } "
            "&& tar -xzf frodock.tgz --strip-components=1 "
            "&& rm -f frodock.tgz",
            targetName=f"{FRODOCK_DIC['name']}_installed")

        installer.addPackage(
            env, dependencies=['wget', 'tar', 'gzip'], default=default)

    @classmethod
    def _addCondaEnvPackage(cls, env, dic, pythonVersion,
                            installCmd=None, default=True):
        activation = cls.getEnvActivationCommand(dic)
        home = cls.getVar(dic['home'])

        installer = InstallHelper(
            dic['name'], packageHome=home, packageVersion=dic['version'])

        installer.getCondaEnvCommand(
            binaryName=dic['name'],
            binaryVersion=dic['version'],
            pythonVersion=pythonVersion)

        if installCmd:
            installer.addCommand(
                f"{activation} && {installCmd}",
                targetName=f"{dic['name']}_installed")

        installer.addCommand(
            f"{activation} && rm -rf {home} && ln -s $CONDA_PREFIX {home}",
            targetName=f"{dic['name']}_symlinked")

        installer.addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addVoromqaPackage(cls, env, default=True):
        cls._addCondaEnvPackage(
            env, VOROMQA_DIC, '3.11',
            f"conda install -y -c bioconda "
            f"voronota={VOROMQA_DIC['version']}",
            default=default)

    @classmethod
    def addFCCPackage(cls, env, default=True):
        installer = InstallHelper(
            FCC_DIC['name'],
            packageHome=cls.getVar(FCC_DIC['home']),
            packageVersion=FCC_DIC['version'])

        installer.getCloneCommand(
            'https://github.com/haddocking/fcc.git',
            binaryFolderName=FCC_DIC['name'],
            targeName=f"{FCC_DIC['name']}_cloned"
        ).addCommand(
            f"cd {FCC_DIC['name']}/src && make",
            targetName=f"{FCC_DIC['name']}_built")

        installer.addPackage(
            env, dependencies=['git', 'make', 'gcc'], default=default)

    @classmethod
    def addProtacModelPackage(cls, env, default=True):
        installer = InstallHelper(
            PROTAC_MODEL_DIC['name'],
            packageHome=cls.getVar(PROTAC_MODEL_DIC['home']),
            packageVersion=PROTAC_MODEL_DIC['version'])

        installer.getCloneCommand(
            'https://github.com/gaoqiweng/PROTAC-Model.git',
            binaryFolderName=PROTAC_MODEL_DIC['name'],
            targeName=f"{PROTAC_MODEL_DIC['name']}_cloned")

        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProtacModelPythonPackage(cls, env, default=True):
        cls._addCondaEnvPackage(
            env, PROTAC_MODEL_PYTHON_DIC, '2.7',
            'conda install -y -c rdkit rdkit=2016.03.3',
            default=default)

    @classmethod
    def addProsettaCPackage(cls, env, default=True):
        installer = InstallHelper(
            PROSETTAC_DIC['name'],
            packageHome=cls.getVar(PROSETTAC_DIC['home']),
            packageVersion=PROSETTAC_DIC['version'])

        installer.getCloneCommand(
            'https://github.com/LondonLab/PRosettaC.git',
            binaryFolderName=PROSETTAC_DIC['name'],
            targeName=f"{PROSETTAC_DIC['name']}_cloned")

        installer.addPackage(env, dependencies=['git'], default=default)

    @classmethod
    def addProsettaCPythonPackage(cls, env, default=True):
        cls._addCondaEnvPackage(
            env, PROSETTAC_PYTHON_DIC, '3.10',
            'conda install -y -c conda-forge rdkit=2023.09.6 '
            'numpy=1.26.4 scikit-learn=1.5.2',
            default=default)

    @classmethod
    def addProsettaCPython2Package(cls, env, default=True):
        cls._addCondaEnvPackage(
            env, PROSETTAC_PYTHON2_DIC, '2.7', default=default)

    # -------------------------------------------------------------------------
    # Installation paths and environment
    # -------------------------------------------------------------------------

    @classmethod
    def _requireToolHome(cls, toolDic, binary=''):
        home = cls.getVar(toolDic['home'])

        if home is None or not os.path.exists(os.path.join(home, binary)):
            raise FileNotFoundError(
                f"{toolDic['home']} does not point to an existing "
                f"{toolDic['name']} installation (got: {home}).")

        return home

    @classmethod
    def _requirePatchdockHome(cls):
        home = cls.getVar(PATCHDOCK_DIC['home'])

        if home is None or not os.path.isdir(home):
            raise FileNotFoundError(
                f"{PATCHDOCK_DIC['home']} is not set. PatchDock must be "
                "installed manually and requires its own registration.")

        return home

    @classmethod
    def _binaryLoads(cls, path):
        if not os.path.exists(path):
            return False

        result = subprocess.run(
            ['ldd', path], capture_output=True, text=True)

        return 'not found' not in result.stdout.lower()

    @classmethod
    def _forceSymlink(cls, source, link):
        pwutils.makePath(os.path.dirname(link))
        pwutils.cleanPath(link)
        pwutils.createAbsLink(source, link)
        return link

    @classmethod
    def prepareFrodockBinDir(cls, targetDir):
        home = cls._requireToolHome(FRODOCK_DIC)
        binDir = os.path.join(targetDir, 'bin')
        os.makedirs(binDir, exist_ok=True)

        for name in FRODOCK_BINARIES:
            intel = os.path.join(home, 'bin', name)
            gcc = os.path.join(home, 'bin', f'{name}_gcc')

            if cls._binaryLoads(intel):
                chosen = intel
            elif cls._binaryLoads(gcc):
                chosen = gcc
            else:
                raise FileNotFoundError(
                    f"Neither {intel} nor {gcc} is usable (checked with ldd).")

            cls._forceSymlink(chosen, os.path.join(binDir, name))

        soapSrc = os.path.join(home, 'bin', 'soap.bin')
        if not os.path.exists(soapSrc):
            raise FileNotFoundError(
                f"{soapSrc} not found under FRODOCK_HOME/bin.")

        cls._forceSymlink(soapSrc, os.path.join(binDir, 'soap.bin'))
        return targetDir

    @classmethod
    def requireRosettaScriptsBinary(cls):
        home = RosettaPlugin.getVar(ROSETTA_DIC['home'])

        if home is None or not os.path.isdir(home):
            raise FileNotFoundError(
                f"{ROSETTA_DIC['home']} does not point to an existing "
                f"Rosetta installation (got: {home}).")

        binDir = os.path.join(home, 'main', 'source', 'bin')

        for build in ROSETTA_SCRIPTS_BUILDS:
            path = os.path.join(
                binDir, f'rosetta_scripts.{build}.linuxgccrelease')
            if os.path.exists(path):
                return path

        others = sorted(glob.glob(
            os.path.join(binDir, 'rosetta_scripts.*')))

        if others:
            return others[0]

        raise FileNotFoundError(
            f"No rosetta_scripts binary found in {binDir}.")

    @classmethod
    def prepareRosettaScriptsShim(cls, targetDir):
        realHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])

        cls._forceSymlink(
            cls.requireRosettaScriptsBinary(),
            os.path.join(
                targetDir, 'main', 'source', 'bin',
                'rosetta_scripts.default.linuxgccrelease'))

        scripts = os.path.join('main', 'source', 'scripts')
        database = os.path.join('main', 'database')

        for relPath, candidates in (
                (scripts, [scripts]),
                (database, [database]),
                ('tools', ['tools', os.path.join('main', 'tools')])):

            sources = [
                os.path.join(realHome, candidate)
                for candidate in candidates
                if os.path.isdir(os.path.join(realHome, candidate))
            ]

            if not sources:
                raise FileNotFoundError(
                    f"{relPath} not found under ROSETTA_HOME ({realHome}).")

            cls._forceSymlink(
                sources[0], os.path.join(targetDir, relPath))

        return targetDir

    @classmethod
    def requireObabelBinary(cls):
        obabelBin = pwchemPlugin.getEnvPath(
            pwchemConstants.OPENBABEL_DIC, innerPath='bin/obabel')

        if not os.path.exists(obabelBin):
            raise FileNotFoundError(
                f"obabel not found at {obabelBin}. Is pwchem's "
                "OpenBabel environment installed?")

        return obabelBin

    @classmethod
    def preparePRosettaCBabelShim(cls, targetDir):
        obabelBin = cls.requireObabelBinary()
        os.makedirs(targetDir, exist_ok=True)

        wrapperPath = os.path.join(targetDir, 'babel')
        activate = pwchemPlugin.getEnvActivationCommand(
            pwchemConstants.OPENBABEL_DIC)

        script = (
            '#!/usr/bin/env bash\n'
            'set -e\n'
            f'{activate}\n'
            f'OBABEL="{os.path.abspath(obabelBin)}"\n'
            'IN="$1"; OUT="$2"; shift 2 || true\n'
            'if [ "$IN" = "$OUT" ]; then\n'
            '  TMP="$(mktemp --suffix=".${OUT##*.}")"\n'
            '  "$OBABEL" "$IN" -O "$TMP" "$@"\n'
            '  mv "$TMP" "$OUT"\n'
            'else\n'
            '  "$OBABEL" "$IN" -O "$OUT" "$@"\n'
            'fi\n'
        )

        with open(wrapperPath, 'w') as f:
            f.write(script)

        os.chmod(wrapperPath, 0o755)
        return targetDir

    @classmethod
    def getProtacModelScript(cls, scriptName=''):
        home = os.path.join(
            cls._requireToolHome(PROTAC_MODEL_DIC),
            PROTAC_MODEL_DIC['name'])

        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getProtacModelPython(cls):
        home = cls._requireToolHome(PROTAC_MODEL_PYTHON_DIC)
        path = os.path.join(home, 'bin', 'python2')

        if not os.path.exists(path):
            raise FileNotFoundError(
                f"python2 not found under PROTAC_MODEL_PYTHON_HOME/bin ({home}).")

        return path

    @classmethod
    def validateInstallation(cls):
        return []

    @classmethod
    def checkProtacModelBinaries(cls):
        frodock = os.path.join(
            cls._requireToolHome(FRODOCK_DIC), 'bin')
        adfr = os.path.join(
            cls._requireToolHome(ADFRSUITE_DIC), 'bin')
        vinaHome = cls._requireToolHome(VINA_DIC)
        voromqa = cls._requireToolHome(VOROMQA_DIC)
        fcc = os.path.join(
            cls._requireToolHome(FCC_DIC), FCC_DIC['name'])

        required = [
            os.path.join(adfr, name)
            for name in (
                'obabel', 'obenergy', 'prepare_ligand',
                'prepare_receptor', 'reduce')
        ]

        required += [
            os.path.join(vinaHome, 'bin', 'vina'),
            os.path.join(voromqa, 'bin', 'voronota-voromqa'),
            os.path.join(fcc, 'src', 'contact_fcc'),
            os.path.join(fcc, 'scripts', 'make_contacts.py'),
            os.path.join(frodock, 'soap.bin'),
            cls.getProtacModelScript(os.path.join('utils', 'frodock.py')),
        ]

        errors = [
            f"{path} not found."
            for path in required if not os.path.exists(path)
        ]

        for name in FRODOCK_BINARIES:
            builds = [
                os.path.join(frodock, name),
                os.path.join(frodock, f'{name}_gcc')
            ]

            if not any(os.path.exists(path) for path in builds):
                errors.append(
                    f"Neither {builds[0]} nor {builds[1]} found.")

        return errors

    @classmethod
    def getProtacModelEnviron(cls, frodockHome=None, needRosetta=False):
        # These names are the environment variables expected by PROTAC-Model.
        # ADFRSuite and Vina paths come from the existing autodock installation.
        environ = {
            'FRODOCK': frodockHome or cls._requireToolHome(FRODOCK_DIC),
            'ADFRSUITE': cls._requireToolHome(ADFRSUITE_DIC),
            'VINA': cls._requireToolHome(VINA_DIC),
            'VOROMQA': cls._requireToolHome(VOROMQA_DIC),
            'FCC': os.path.join(
                cls._requireToolHome(FCC_DIC), FCC_DIC['name']),
            'PROTAC_MODEL_HOME': cls.getProtacModelScript(),
        }

        if needRosetta:
            rosettaHome = RosettaPlugin.getVar(ROSETTA_DIC['home'])

            if rosettaHome is None or not os.path.isdir(rosettaHome):
                raise FileNotFoundError(
                    f"{ROSETTA_DIC['home']} is not set or does not exist "
                    f"(got: {rosettaHome}).")

            environ['ROSETTA'] = rosettaHome

        return {
            name: os.path.abspath(path)
            for name, path in environ.items()
        }

    @classmethod
    def getPRosettaCScript(cls, scriptName=''):
        home = os.path.join(
            cls._requireToolHome(PROSETTAC_DIC), PROSETTAC_DIC['name'])

        return os.path.join(home, scriptName) if scriptName else home

    @classmethod
    def getPRosettaCEnviron(cls, rosettaHome, obDir):
        python2Home = cls._requireToolHome(PROSETTAC_PYTHON2_DIC)

        return {
            'PATCHDOCK': cls._requirePatchdockHome(),
            'OB': os.path.abspath(obDir),
            'SCRIPTS_FOL': os.path.abspath(
                cls.getPRosettaCScript()) + os.sep,
            'ROSETTA_FOL': os.path.abspath(rosettaHome),
            'PATH': os.environ.get('PATH', '') + os.pathsep
                    + os.path.join(python2Home, 'bin'),
        }

    @classmethod
    def getPluginScript(cls, scriptName):
        return os.path.join(_PLUGIN_DIR, 'scripts', scriptName)

    # -------------------------------------------------------------------------
    # Running external scripts
    # -------------------------------------------------------------------------

    @classmethod
    def runCondaScript(cls, protocol, scriptPath, args, condaDic,
                       extraEnvDict=None, cwd=None):
        program = (
            f'{cls.getEnvActivationCommand(condaDic)} && '
            f'python "{os.path.abspath(scriptPath)}"'
        )

        cls.runProgram(
            protocol, program, args,
            extraEnvDict=extraEnvDict, cwd=cwd)

    @classmethod
    def getEnviron(cls):
        return pwutils.Environ(os.environ)

    @classmethod
    def runProgram(cls, protocol, program, args=None,
                   extraEnvDict=None, cwd=None):
        env = cls.getEnviron()

        if extraEnvDict is not None:
            env.update(extraEnvDict)

        protocol.runJob(
            program, args, env=env,
            cwd=os.path.abspath(cwd) if cwd else cwd,
            numberOfMpi=1)