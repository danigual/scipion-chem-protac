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

"""
Wrapper around PROTAC-Model (Gao et al., Zhejiang University,
https://github.com/gaoqiweng/PROTAC-Model), an integrative modeling pipeline for
PROTAC-mediated protein-protein ternary complexes:

  1. FRODOCK: global rigid-body protein-protein docking, guided by a site point on the
     receptor interface.
  2. Filtering of the FRODOCK poses by compatibility with the PROTAC SMILES (and,
     optionally, up to two E3 ligand SDF files, for PROTACs with more than one possible
     anchoring orientation).
  3. Optional RosettaDock refinement of the filtered poses (slower, more accurate).

The pipeline logic is not reimplemented: each step launches the Python 2 driver
run_protac_model.py, which calls PROTAC-Model's own functions.
"""

import json, csv
import os

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message, cleanPath
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SetOfAtomStructsChem
from pwchem.utils import cleanPDB, convertToSdf
from pwem.convert import cifToPdb

from protac import Plugin, PROTAC_MODEL_PYTHON_DIC
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC
from protac.utils.molecules import findMolByName, listMolNames


class ProtPROTACModel(EMProtocol):
    """
    Models a PROTAC-mediated protein-protein ternary complex (receptor + target,
    bridged by a PROTAC) using the PROTAC-Model pipeline (Gao et al.): FRODOCK global
    docking guided by a receptor-interface site point, filtering of the resulting poses
    by compatibility with the PROTAC linker geometry, and optional RosettaDock
    refinement of the surviving poses.
    """
    _label = 'PROTAC ternary complex modeling (PROTAC-Model)'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Proteins')
        group.addParam('inputReceptor', params.PointerParam, pointerClass='AtomStruct',
                       label='Receptor structure',
                       help='Larger of the two proteins. Passed '
                            'as the FRODOCK receptor. Per PROTAC-Model\'s own requirements, '
                            'this structure should already include its bound small-molecule '
                            'warhead (as HETATM records). Unrelated heteroatoms are stripped '
                            'automatically when its warhead is named below.')
        #todo see how i do this, for the site thing maybe do it internally directly
        group.addParam('receptorLigandName', params.StringParam,
                       label='Receptor warhead residue name',
                       help='Residue name of the bound warhead in the receptor, as it appears '
                            'in the PDB HETATM records (e.g. "B96"). Every other heteroatom '
                            'is then removed.\n\n'
                            'PROTAC-Model reads every HETATM record as part of the warhead, '
                            'so a leftover sulfate or glycerol would shift the anchoring '
                            'point and give wrong models with normal-looking scores. Leave '
                            'empty only if the warhead is already the only heteroatom.')
        group.addParam('inputTarget', params.PointerParam, pointerClass='AtomStruct',
                       label='Target structure',
                       help='Smaller of the two proteins (typically the protein of '
                            'interest, POI). Passed as the FRODOCK docking target. Per '
                            'PROTAC-Model\'s own requirements, this structure should '
                            'already include its bound small-molecule warhead (as HETATM '
                            'records). Unrelated heteroatoms are stripped automatically '
                            'when its warhead is named below.')
        group.addParam('targetLigandName', params.StringParam,
                       label='Target warhead residue name',
                       help='Residue name of the bound warhead in the target, as it '
                            'appears in the PDB HETATM records. Same purpose and same '
                            'caveats as "Receptor warhead residue name" above, for the '
                            'other protein.')
        #todo do internally
        group.addParam('siteCoords', params.StringParam,
                       label='Receptor interface site (X,Y,Z)',
                       help='Point on the receptor where the ternary complex is expected '
                            'to form; FRODOCK restricts its search around it. The wizard '
                            'fills it from the receptor warhead named above.\n\n'
                            "It must be the centroid of the receptor's warhead HETATM "
                            'records, measured on the cleaned structure (warhead only) '
                            'when "Receptor warhead residue name" is set.')

        group = form.addGroup('PROTAC')
        group.addParam('protacMolecule', params.PointerParam,
                       pointerClass='SetOfSmallMolecules',
                       label='PROTAC molecule: ',
                       help='Full PROTAC molecule (E3 ligand - linker - '
                            'POI warhead), used to filter FRODOCK poses by compatibility '
                            'with the linker geometry.')

        group.addParam('receptorLigandMol', params.PointerParam,
                       pointerClass='SetOfSmallMolecules',
                       label='Receptor ligand: ',
                       help='Docked receptor ligand. ')

        group.addParam('targetLigandMol', params.PointerParam,
                       pointerClass='SetOfSmallMolecules',
                       label='Target ligand: ',
                       help='Docked target ligand. ')

        group = form.addGroup('Docking parameters')
        group.addParam('doRefine', params.BooleanParam, default=False,
                       label='Refine with RosettaDock',
                       help='Refine the filtered FRODOCK poses with RosettaDock. Improves '
                            'accuracy but is considerably slower than FRODOCK + filtering '
                            'alone.\n\n'
                            'Requires more than 10 threads ("Threads" below): '
                            "PROTAC-Model's own refinement code crashes outright with 10 "
                            'or fewer.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllStepsOld(self):
        # Scipion compares each step's arguments on its own when continuing a run, so every
        # step also gets a key with its predecessors' arguments: changing an input then
        # reruns the step that uses it and everything after it.
        convertArgs = [self.inputReceptor.get().getFileName(),
                       self.receptorLigandName.get() or '',
                       self.inputTarget.get().getFileName(),
                       self.targetLigandName.get() or '',
                       *self._getE3LigandArgs()]
        e3lig1, e3lig2 = self._getE3LigandFiles()
        dockArgs = [self._getSiteArg(), (self.protacSmiles.get() or '').strip(), e3lig1, e3lig2]
        # Threads only change the results when refining (they set Rosetta's nstruct).
        filterArgs = [self._getLigLocateNum(), self._getSmiArg(self.targetLigandSmiles),
                      self._getSmiArg(self.receptorLigandSmiles)]

        convertId = self._insertFunctionStep(self.convertInputStep, *convertArgs,
                                             prerequisites=[])
        dockId = self._insertFunctionStep(self.frodockStep, *dockArgs,
                                          self._upstreamKey(convertArgs),
                                          prerequisites=[convertId])
        lastId = self._insertFunctionStep(self.filterPosesStep, *filterArgs,
                                          self._upstreamKey(convertArgs, dockArgs),
                                          prerequisites=[dockId])
        refineArgs = []
        if self.doRefine.get():
            refineArgs = [self.numberOfThreads.get()] + filterArgs
            lastId = self._insertFunctionStep(
                self.refineStep, *refineArgs,
                self._upstreamKey(convertArgs, dockArgs, filterArgs),
                prerequisites=[lastId])

        self._insertFunctionStep(
            self.createOutputStep,
            self._upstreamKey(convertArgs, dockArgs, filterArgs, refineArgs,
                              [self.doRefine.get()]),
            prerequisites=[lastId])

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.frodockStep)
        self._insertFunctionStep(self.filterPosesStep)

    def convertInputStep(self):
        for structure, basename in (
                (self.inputReceptor.get().getFileName(), 'receptor'),
                (self.inputTarget.get().getFileName(), 'target')
        ):
            file = structure

            if file.endswith('.cif'):
                inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
                cifToPdb(file, inpFile)
            else:
                inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
                shutil.copy2(file, inpFile)

        ligands = [
            ('receptor', self.receptorLigandMol.get().getFirstItem()),
            ('target', self.targetLigandMol.get().getFirstItem()),
            ('protac', self.protacMolecule.get().getFirstItem())
        ]

        extraDir = os.path.abspath(self._getExtraPath())
        smilesList = os.path.join(extraDir, 'ligands.txt')

        with open(smilesList, 'w') as f:
            for name, mol in ligands:
                molFile = (
                    mol.getPoseFile()
                    if mol.getPoseFile()
                    else mol.getFileName()
                )
                molFile = os.path.abspath(molFile)

                # Convert MOL files to SDF before RDKit processing
                if os.path.splitext(molFile)[1].lower() == '.mol':
                    sdfFile = os.path.join(extraDir, f'{name}_ligand.sdf')

                    pwchemPlugin.runScript(
                        self,
                        'obabel_IO.py',
                        [
                            '-i', molFile,
                            '-of', 'sdf',
                            '-o', f'{name}_ligand',
                            '-od', extraDir
                        ],
                        env=OPENBABEL_DIC,
                        cwd=self._getPath()
                    )

                    molFile = sdfFile

                f.write(f'{molFile}\n')

        # Generate SMILES CSV using RDKit
        smilesCsv = os.path.join(extraDir, 'ligand_smiles.csv')

        pwchemPlugin.runScript(
            self,
            'rdkit_IO.py',
            [
                '-i', smilesList,
                '-of', 'smiles_csv',
                '-o', 'ligand_smiles',
                '-od', extraDir
            ],
            env=RDKIT_DIC,
            cwd=self._getPath()
        )

        with open(smilesCsv, newline='') as f:
            reader = csv.DictReader(f)
            smilesColumn = next(
                (
                    col for col in (reader.fieldnames or [])
                    if 'smiles' in col.lower()
                ),
                None
            )

            if smilesColumn is None:
                raise ValueError(
                    f'No SMILES column found in {smilesCsv}. '
                    f'CSV columns: {reader.fieldnames}'
                )

            rows = list(reader)

        self.receptorLigandSmiles = rows[0][smilesColumn].strip()
        self.targetLigandSmiles = rows[1][smilesColumn].strip()
        self.protacSmiles = rows[2][smilesColumn].strip()

        print(f'Receptor ligand SMILES: {self.receptorLigandSmiles}')
        print(f'Target ligand SMILES: {self.targetLigandSmiles}')
        print(f'PROTAC SMILES: {self.protacSmiles}')

        #todo see if i can exclude this and call it directly instead of copying
        Plugin.prepareFrodockBinDir(self._getFrodockBinDir())

    def frodockStep(self):
        frodockDir = self._getExtraPath('frodock')
        os.makedirs(frodockDir, exist_ok=True)

        receptorFile = os.path.abspath(self._getReceptorFile())
        targetFile = os.path.abspath(self._getTargetFile())

        args = (
            f'--phase frodock '
            f'--receptor "{receptorFile}" '
            f'--target "{targetFile}" '
            f'--smiles "{self.protacSmiles}" '
            f'--site={self._getSiteArg()}'
        )

        Plugin.runCondaScript(
            self,
            Plugin.getPluginScript('run_protac_model.py'),
            args,
            PROTAC_MODEL_PYTHON_DIC,
            extraEnvDict=Plugin.getProtacModelEnviron(),
            cwd=frodockDir
        )

    def filterPosesStep(self):
        ligLocateNum = 1
        cleanPath(self._getExtraPath('frodock_results'))

        args = (f'--phase filter --cpu {self.numberOfThreads.get()} '
               f'--lig-locate-num {ligLocateNum} '
               f'--target-smi "{self.targetLigandSmiles}" --rec-smi "{self.receptorLigandSmiles}"')
        Plugin.runCondaScript(self, Plugin.getPluginScript('run_protac_model.py'), args,
                              PROTAC_MODEL_PYTHON_DIC,
                              extraEnvDict=Plugin.getProtacModelEnviron(
                                  frodockHome=self._getFrodockBinDir()),
                              cwd=self._getExtraPath('frodock'))

    def refineStep(self, cpu, ligLocateNum, targetSmi, recSmi, upstream):
        """ Optional RosettaDock refinement (ros.rosetta). """
        # Next to extra/frodock/: ros.rosetta() reaches it through '../frodock/'.
        rosettaDir = self._getExtraPath('rosetta')
        os.makedirs(rosettaDir, exist_ok=True)
        # Same reason as in filterPosesStep.
        cleanPath(self._getExtraPath('rosetta_results'))

        args = (f'--phase refine --cpu {cpu} --lig-locate-num {ligLocateNum} '
               f'--target-smi "{targetSmi}" --rec-smi "{recSmi}"')
        Plugin.runCondaScript(self, Plugin.getPluginScript('run_protac_model.py'), args,
                              PROTAC_MODEL_PYTHON_DIC,
                              extraEnvDict=Plugin.getProtacModelEnviron(
                                  frodockHome=self._getFrodockBinDir(), needRosetta=True),
                              cwd=rosettaDir)

    def createOutputStep(self, upstream):
        """ One AtomStruct per final model, with its VoroMQA score (lower is better) as
        '_score'. """
        if self.doRefine.get():
            resultsDir = self._getExtraPath('rosetta_results', 'all')
            resultsFile = os.path.join(resultsDir, 'results_rosetta.txt')
        else:
            resultsDir = self._getExtraPath('frodock_results', 'all')
            resultsFile = os.path.join(resultsDir, 'results_frodock.txt')

        # 0 poses must fail the protocol instead of finishing green with an empty set.
        if not os.path.exists(resultsFile):
            raise RuntimeError(
                f'{resultsFile} was not produced. This usually means PROTAC-Model found '
                '0 compatible poses: either the receptor/target warheads are not exact '
                'substructures of the PROTAC SMILES (PROTAC-Model needs '
                'GetSubstructMatch to succeed for both), or "Receptor interface site '
                '(X,Y,Z)" is not the centroid of the receptor\'s warhead HETATM records.')

        outputSet = SetOfAtomStructsChem().create(self._getPath())
        with open(resultsFile) as f:
            for line in f:
                # "<pose_id> <score>", no header.
                fields = line.split()
                if len(fields) < 2:
                    continue
                poseId, score = fields[:2]
                pdbFile = os.path.join(resultsDir, f'model_merge_{poseId}.pdb')
                if not os.path.exists(pdbFile):
                    self.info(f'Skipping pose {poseId}: {pdbFile} not found.')
                    continue
                atomStruct = AtomStruct(filename=pdbFile)
                atomStruct.setObjLabel(f'model_merge_{poseId}')
                # AtomStruct has no score field; '_score' follows pwchem's convention.
                atomStruct._score = pwobj.Float(float(score))
                outputSet.append(atomStruct)

        if len(outputSet) == 0:
            raise RuntimeError(
                f'{resultsFile} exists but lists no usable pose (either it is empty, or '
                'every model_merge_*.pdb it references is missing). This usually means '
                'PROTAC-Model found 0 compatible poses: either the receptor/target '
                'warheads are not exact substructures of the PROTAC SMILES, or '
                '"Receptor interface site (X,Y,Z)" is not the centroid of the '
                "receptor's warhead HETATM records.")

        self._defineOutputs(outputTernaryModels=outputSet)
        self._defineSourceRelation(self.inputReceptor, outputSet)
        self._defineSourceRelation(self.inputTarget, outputSet)
        # The E3 conformers also decide which poses survive, so they go in the provenance.
        if self._useE3Ligands():
            self._defineSourceRelation(self.e3Ligands, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        summary = []
        if self.isFinished() and hasattr(self, 'outputTernaryModels'):
            models = self.outputTernaryModels
            if len(models):
                bestScore = min(model._score.get() for model in models)
                summary.append(f'Generated {len(models)} ternary complex model(s); best '
                               f'(most negative) interface score: {bestScore:.2f}.')
        return summary

    def _citations(self):
        return []

    # --------------------------- UTILS functions ------------------------------
    def getE3LigandNames(self):
        """ Molecule names in the E3 conformer set, listed by the name wizards. """
        return listMolNames(self.e3Ligands.get())

    @staticmethod
    def _cleanReceptorOrTarget(inFile, outFile, ligandName):
        """ Water-stripped copy of inFile keeping only the named warhead, if any, among
        heteroatoms. Upper-cased because PDB residue names are. """
        if ligandName and ligandName.strip():
            cleanPDB(inFile, outFile, waters=True, hetatm=True,
                    het2keep=[ligandName.strip().upper()])
        else:
            cleanPDB(inFile, outFile, waters=True, hetatm=False)

    def _getSiteCoords(self):
        """ Parses siteCoords ("X,Y,Z") into a tuple of 3 floats, or None if malformed. """
        parts = (self.siteCoords.get() or '').strip().split(',')
        if len(parts) != 3:
            return None
        try:
            return tuple(float(p.strip()) for p in parts)
        except ValueError:
            return None

    def _getSiteArg(self):
        """ siteCoords rebuilt from the parsed values, so stray spaces can't split the
        shell argument. Raw text if malformed (_validate reports it). """
        coords = self._getSiteCoords()
        if coords is None:
            return self.siteCoords.get() or ''
        return '%.4f,%.4f,%.4f' % coords

    @staticmethod
    def _upstreamKey(*argLists):
        return json.dumps([arg for args in argLists for arg in args])

    def _useE3Ligands(self):
        return (self.e3Ligands.get() is not None and bool(self.e3Ligand1Name.get())
                and bool(self.e3Ligand2Name.get()))

    def _getE3LigandArgs(self):
        """ Set id and both conformer names, or empty values without conformers. """
        if not self._useE3Ligands():
            return [0, '', '']
        return [self.e3Ligands.get().getObjId(), self.e3Ligand1Name.get().strip(),
                self.e3Ligand2Name.get().strip()]

    def _getE3LigandFiles(self):
        """ The two E3 conformer SDFs written by convertInputStep, or two empty strings:
        PROTAC-Model uses them only as a pair. """
        if not self._useE3Ligands():
            return '', ''
        return self._getExtraPath('e3_ligand_1.sdf'), self._getExtraPath('e3_ligand_2.sdf')

    def _getReceptorFile(self):
        """ Recomputed rather than cached: each step may run in its own process. """
        return self._getExtraPath('receptor.pdb')

    def _getTargetFile(self):
        """ See _getReceptorFile. """
        return self._getExtraPath('target.pdb')

    def _getFrodockBinDir(self):
        """ FRODOCK binary shim built by convertInputStep. """
        return self._getExtraPath('frodock_bin')

    def _getLigLocateNum(self):
        """ PROTAC-Model's lig_locate_num: 2 with both E3 conformers, 1 otherwise. """
        return 2 if all(self._getE3LigandFiles()) else 1

    @staticmethod
    def _getSmiArg(smilesParam):
        """ 'none' is how PROTAC-Model means "no SMILES given". """
        value = (smilesParam.get() or '').strip()
        return value or 'none'


