# **************************************************************************
# *
# * Authors:     Daniel Gutiérrez (Internship student, github user: danigual)
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
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

# Scipion em imports
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.protocol import getProtocolFromDb
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules, ProtExtractLigands, ConvertStructures
from pwchem.utils import assertHandle

# Plugin imports
from .. import Plugin
from ..protocols import ProtPROTACModel

# example_1 of the PROTAC-Model repository (VHL + SMARCA2 bridged by its PROTAC). It
# ships with the PROTAC-Model checkout installed by this plugin.
EXAMPLE_DIR = 'example_1'
RECEPTOR_FILE, TARGET_FILE = 'receptor.pdb', 'target.pdb'
PROTAC_PATTERN = 'protac.smi'
SITE_FILE = 'site_info.txt'

# Residue name of the warhead bound to both proteins in example_1.
WARHEAD_NAME = 'FWZ'

SDF_FORMAT, OBABEL_MANAGER = 2, 1

# FRODOCK docking plus pose filtering takes around 2 h, and the optional refinement
# considerably longer. _waitOutput gives up after timeOut seconds, which also has to
# cover the time the protocol spends waiting in the queue, if one is used.
SLEEP_TIME = 60
MODEL_TIMEOUT, REFINE_TIMEOUT = 4 * 3600, 12 * 3600


class TestPROTACModel(BaseTest):
    """ PROTAC-Model ternary complex modeling: FRODOCK docking of the two proteins
    around the receptor interface site, plus filtering of the poses by compatibility
    with the PROTAC. Run on example_1 of the PROTAC-Model repository. """

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.exampleDir = cls._getExampleDir()
        if cls.exampleDir is None:
            return

        # Imports are launched together, then each stage waits for the outputs it needs
        cls._runImportProteins()
        cls._runImportPROTAC()
        cls._waitOutput(cls.protImportReceptor, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportTarget, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportPROTAC, 'outputSmallMolecules', sleepTime=5)

        cls._runExtractWarheads()
        cls._waitOutput(cls.protExtractRec, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protExtractTarget, 'outputSmallMolecules', sleepTime=5)

        cls._runConvertWarheads()
        cls._waitOutput(cls.protConvertRec, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protConvertTarget, 'outputSmallMolecules', sleepTime=5)

    # --------------------------- Input protocols ------------------------------
    @classmethod
    def _runImportProteins(cls):
        cls.protImportReceptor = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.getExampleFile(RECEPTOR_FILE))
        cls.proj.launchProtocol(cls.protImportReceptor, wait=False)

        cls.protImportTarget = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.getExampleFile(TARGET_FILE))
        cls.proj.launchProtocol(cls.protImportTarget, wait=False)

    @classmethod
    def _runImportPROTAC(cls):
        # make3d writes the molecule as SDF, the format the protocol reads.
        cls.protImportPROTAC = cls.newProtocol(
            ProtChemImportSmallMolecules,
            singleFiles=True, filesPath=cls.exampleDir, filesPattern=PROTAC_PATTERN,
            useManager=OBABEL_MANAGER, make3d=True)

        cls.proj.launchProtocol(cls.protImportPROTAC, wait=False)

    @classmethod
    def _runExtractWarheads(cls):
        cls.protExtractRec = cls._runExtractLigand(cls.protImportReceptor)
        cls.protExtractTarget = cls._runExtractLigand(cls.protImportTarget)

    @classmethod
    def _runExtractLigand(cls, inputProt):
        protExtract = cls.newProtocol(ProtExtractLigands)

        protExtract.inputStructure.set(inputProt)
        protExtract.inputStructure.setExtended('outputPdb')

        cls.proj.launchProtocol(protExtract, wait=False)
        return protExtract

    @classmethod
    def _runConvertWarheads(cls):
        cls.protConvertRec = cls._runConvertToSdf(cls.protExtractRec)
        cls.protConvertTarget = cls._runConvertToSdf(cls.protExtractTarget)

    @classmethod
    def _runConvertToSdf(cls, inputProt):
        # The extracted warheads are PDB pose files; the protocol only reads SDF or MOL.
        protConvert = cls.newProtocol(
            ConvertStructures,
            outputFormatSmall=SDF_FORMAT, usePose=True, useManager=OBABEL_MANAGER)

        protConvert.inputObject.set(inputProt)
        protConvert.inputObject.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(protConvert, wait=False)
        return protConvert

    # --------------------------- Tested protocol ------------------------------
    def _runPROTACModel(self, doRefine=False, threads=4):
        protModel = self.newProtocol(
            ProtPROTACModel,
            receptorLigandName=WARHEAD_NAME, targetLigandName=WARHEAD_NAME,
            siteCoords=self.getSiteCoords(),
            doRefine=doRefine, numberOfThreads=threads)

        protModel.inputReceptor.set(self.protImportReceptor)
        protModel.inputReceptor.setExtended('outputPdb')

        protModel.inputTarget.set(self.protImportTarget)
        protModel.inputTarget.setExtended('outputPdb')

        protModel.receptorLigandMol.set(self.protConvertRec)
        protModel.receptorLigandMol.setExtended('outputSmallMolecules')

        protModel.targetLigandMol.set(self.protConvertTarget)
        protModel.targetLigandMol.setExtended('outputSmallMolecules')

        protModel.protacMolecule.set(self.protImportPROTAC)
        protModel.protacMolecule.setExtended('outputSmallMolecules')

        self.proj.launchProtocol(protModel, wait=False)
        return protModel

    def test(self):
        self._checkExampleDir()

        protModel = self._runPROTACModel()
        self._waitOutput(protModel, 'outputTernaryModels',
                         sleepTime=SLEEP_TIME, timeOut=MODEL_TIMEOUT)
        protModel = self._updateProt(protModel)

        models = getattr(protModel, 'outputTernaryModels', None)
        assertHandle(self.assertIsNotNone, models, cwd=protModel.getWorkingDir(),
                     message='No ternary complex models were generated')
        assertHandle(self.assertGreater, len(models), 0, cwd=protModel.getWorkingDir(),
                     message='The output set of ternary complex models is empty')

    # --------------------------- UTILS functions ------------------------------
    @classmethod
    def _updateProt(cls, prot):
        """ Reload a protocol from its own database. _waitOutput only reloads its local
        copy, so the protocol object held here still lacks the outputs it registered. """
        return getProtocolFromDb(prot.getProject().path, prot.getDbPath(), prot.getObjId())

    @classmethod
    def _getExampleDir(cls):
        """ example_1 inside the installed PROTAC-Model, or None when PROTAC-Model is
        not installed. """
        try:
            exampleDir = Plugin.getProtacModelScript(EXAMPLE_DIR)
        except FileNotFoundError:
            return None

        return exampleDir if os.path.isdir(exampleDir) else None

    def _checkExampleDir(self):
        if self.exampleDir is None:
            self.skipTest('PROTAC-Model is not installed, so its example data '
                          '(example_1) is not available')

    @classmethod
    def getExampleFile(cls, fileName):
        return os.path.join(cls.exampleDir, fileName)

    @classmethod
    def getSiteCoords(cls):
        """ Receptor interface site shipped with the example, as "X,Y,Z". """
        with open(cls.getExampleFile(SITE_FILE)) as f:
            return f.read().strip()
