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

import os

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message, cleanPath
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SetOfAtomStructsChem
from pwchem.utils import cleanPDB

from protac import Plugin, PROTAC_MODEL_PYTHON_DIC


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
                       label='Receptor structure', allowsNull=False,
                       help='Larger of the two proteins. Passed '
                            'as the FRODOCK receptor. Per PROTAC-Model\'s own requirements, '
                            'this structure should already include its bound small-molecule '
                            'warhead (as HETATM records). Unrelated heteroatoms are stripped '
                            'automatically when its warhead is named below.')
        group.addParam('receptorLigandName', params.StringParam, allowsNull=True,
                       label='Receptor warhead residue name',
                       help='Residue name of the bound warhead in the receptor, as it appears '
                            'in the PDB HETATM records (e.g. "B96"). Every other heteroatom '
                            'is then removed.\n\n'
                            'PROTAC-Model reads every HETATM record as part of the warhead, '
                            'so a leftover sulfate or glycerol would shift the anchoring '
                            'point and give wrong models with normal-looking scores. Leave '
                            'empty only if the warhead is already the only heteroatom.')
        group.addParam('inputTarget', params.PointerParam, pointerClass='AtomStruct',
                       label='Target structure', allowsNull=False,
                       help='Smaller of the two proteins (typically the protein of '
                            'interest, POI). Passed as the FRODOCK docking target. Per '
                            'PROTAC-Model\'s own requirements, this structure should '
                            'already include its bound small-molecule warhead (as HETATM '
                            'records). Unrelated heteroatoms are stripped automatically '
                            'when its warhead is named below.')
        group.addParam('targetLigandName', params.StringParam, allowsNull=True,
                       label='Target warhead residue name',
                       help='Residue name of the bound warhead in the target, as it '
                            'appears in the PDB HETATM records. Same purpose and same '
                            'caveats as "Receptor warhead residue name" above, for the '
                            'other protein.')
        group.addParam('siteCoords', params.StringParam, allowsNull=False,
                       label='Receptor interface site (X,Y,Z)',
                       help='Point on the receptor where the ternary complex is expected '
                            'to form; FRODOCK restricts its search around it.\n\n'
                            "It must be the centroid of the receptor's warhead HETATM "
                            'records, measured on the cleaned structure (warhead only) '
                            'when "Receptor warhead residue name" is set.')

        group = form.addGroup('PROTAC')
        group.addParam('protacSmiles', params.StringParam, allowsNull=False,
                       label='PROTAC SMILES',
                       help='SMILES of the full PROTAC molecule (E3 ligand - linker - '
                            'POI warhead), used to filter FRODOCK poses by compatibility '
                            'with the linker geometry.')
        group.addParam('receptorLigandSmiles', params.StringParam, allowsNull=True,
                       label='Receptor-bound ligand SMILES (optional)',
                       help='SMILES of the small-molecule warhead already bound in the '
                            'receptor structure (equivalent to PROTAC-Model\'s -irsmi). '
                            'RDKit can sometimes fail to assign the correct bonds for this '
                            'ligand when reading it straight from the PDB HETATM records; '
                            'providing its SMILES avoids that failure mode.')
        group.addParam('targetLigandSmiles', params.StringParam, allowsNull=True,
                       label='Target-bound ligand SMILES (optional)',
                       help='SMILES of the small-molecule warhead already bound in the '
                            'target structure (equivalent to PROTAC-Model\'s -itsmi). Same '
                            'purpose as "Receptor-bound ligand SMILES", for the other '
                            'protein.')
        group.addParam('e3Ligand1', params.PointerParam, pointerClass='SmallMolecule',
                       allowsNull=True, label='E3 ligand conformer 1 (optional)',
                       help='SDF of one possible bound conformation/orientation of the '
                            'E3-binding warhead. Optional: only needed to disambiguate '
                            'PROTACs whose E3 ligand can anchor in more than one way '
                            '(e.g. thalidomide-based degraders).')
        group.addParam('e3Ligand2', params.PointerParam, pointerClass='SmallMolecule',
                       allowsNull=True, label='E3 ligand conformer 2 (optional)',
                       condition='e3Ligand1 is not None',
                       help='SDF of a second possible bound conformation/orientation of '
                            'the E3-binding warhead, as an alternative to conformer 1.')

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
    def _insertAllSteps(self):
        convertId = self._insertFunctionStep(self.convertInputStep, prerequisites=[])
        dockId = self._insertFunctionStep(self.frodockStep, prerequisites=[convertId])
        filterId = self._insertFunctionStep(self.filterPosesStep, prerequisites=[dockId])

        lastId = filterId
        if self.doRefine.get():
            lastId = self._insertFunctionStep(self.refineStep, prerequisites=[filterId])

        self._insertFunctionStep(self.createOutputStep, prerequisites=[lastId])

    def convertInputStep(self):
        """ Writes the cleaned receptor/target PDBs: no water, and only the named warhead
        among heteroatoms when its residue name is given. """
        self._cleanReceptorOrTarget(self.inputReceptor.get().getFileName(),
                                    self._getReceptorFile(), self.receptorLigandName)
        self._cleanReceptorOrTarget(self.inputTarget.get().getFileName(),
                                    self._getTargetFile(), self.targetLigandName)

        # Picks the FRODOCK build (intel/gcc) that actually loads on this machine;
        # PROTAC-Model only checks that the intel binary exists.
        Plugin.prepareFrodockBinDir(self._getFrodockBinDir())

    def frodockStep(self):
        """ FRODOCK global docking (fro.frodock). The frodock and filter phases share
        extra/frodock/, since the original works with paths relative to it. """
        frodockDir = self._getExtraPath('frodock')
        os.makedirs(frodockDir, exist_ok=True)

        # Absolute (the driver runs from another cwd) and quoted (it goes through a shell).
        receptorFile = os.path.abspath(self._getReceptorFile())
        targetFile = os.path.abspath(self._getTargetFile())
        # Rebuilt from the parsed values so stray spaces can't split the argument.
        x, y, z = self._getSiteCoords()
        site = f'{x:.4f},{y:.4f},{z:.4f}'
        args = (f'--phase frodock --receptor "{receptorFile}" --target "{targetFile}" '
               f'--smiles "{self.protacSmiles.get().strip()}" --site={site}')

        # Both conformers or neither, matching _getLigLocateNum.
        if self.e3Ligand1.get() is not None and self.e3Ligand2.get() is not None:
            e3lig1 = os.path.abspath(self.e3Ligand1.get().getFileName())
            e3lig2 = os.path.abspath(self.e3Ligand2.get().getFileName())
            args += f' --e3lig1 "{e3lig1}" --e3lig2 "{e3lig2}"'

        Plugin.runCondaScript(Plugin.getPluginScript('run_protac_model.py'), args,
                              PROTAC_MODEL_PYTHON_DIC,
                              extraEnvDict=Plugin.getProtacModelEnviron(
                                  frodockHome=self._getFrodockBinDir()),
                              cwd=frodockDir)

    def filterPosesStep(self):
        """ Filters the FRODOCK poses by compatibility with the PROTAC
        (fro.filter_frodock). """
        # The original never clears old results, and pose ids restart at 1 on every run,
        # so a relaunch could pair a stale model with a fresh score.
        cleanPath(self._getExtraPath('frodock_results'))

        targetSmi = self._getSmiArg(self.targetLigandSmiles)
        recSmi = self._getSmiArg(self.receptorLigandSmiles)
        args = (f'--phase filter --cpu {self.numberOfThreads.get()} '
               f'--lig-locate-num {self._getLigLocateNum()} '
               f'--target-smi "{targetSmi}" --rec-smi "{recSmi}"')

        Plugin.runCondaScript(Plugin.getPluginScript('run_protac_model.py'), args,
                              PROTAC_MODEL_PYTHON_DIC,
                              extraEnvDict=Plugin.getProtacModelEnviron(
                                  frodockHome=self._getFrodockBinDir()),
                              cwd=self._getExtraPath('frodock'))

    def refineStep(self):
        """ Optional RosettaDock refinement (ros.rosetta). """
        # Next to extra/frodock/: ros.rosetta() reaches it through '../frodock/'.
        rosettaDir = self._getExtraPath('rosetta')
        os.makedirs(rosettaDir, exist_ok=True)
        # Same reason as in filterPosesStep.
        cleanPath(self._getExtraPath('rosetta_results'))

        targetSmi = self._getSmiArg(self.targetLigandSmiles)
        recSmi = self._getSmiArg(self.receptorLigandSmiles)
        args = (f'--phase refine --cpu {self.numberOfThreads.get()} '
               f'--lig-locate-num {self._getLigLocateNum()} '
               f'--target-smi "{targetSmi}" --rec-smi "{recSmi}"')

        Plugin.runCondaScript(Plugin.getPluginScript('run_protac_model.py'), args,
                              PROTAC_MODEL_PYTHON_DIC,
                              extraEnvDict=Plugin.getProtacModelEnviron(
                                  frodockHome=self._getFrodockBinDir()),
                              cwd=rosettaDir)

    def createOutputStep(self):
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
                if not line.strip():
                    continue
                poseId, score = line.split()[:2]
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
        if self.e3Ligand1.get() is not None:
            self._defineSourceRelation(self.e3Ligand1, outputSet)
        if self.e3Ligand2.get() is not None:
            self._defineSourceRelation(self.e3Ligand2, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []

        if self._getSiteCoords() is None:
            errors.append('"Receptor interface site (X,Y,Z)" must be 3 comma-separated '
                          f'numbers, e.g. "12.3,-4.5,6.7". Got: "{self.siteCoords.get()}".')

        if self.e3Ligand2.get() is not None and self.e3Ligand1.get() is None:
            errors.append('"E3 ligand conformer 2" was set without "E3 ligand conformer 1". '
                          'Set conformer 1 first, or clear conformer 2.')

        # generate_rosetta_para() only defines gen_conf_num when cpu > 10, so refinement
        # would crash after hours of docking and filtering.
        if self.doRefine.get() and self.numberOfThreads.get() <= 10:
            errors.append('"Refine with RosettaDock" needs more than 10 threads: '
                          'PROTAC-Model\'s own refinement code (generate_rosetta_para) '
                          'crashes with cpu<=10 (UnboundLocalError on gen_conf_num), and '
                          'that number also sets mpirun\'s -np and Rosetta\'s nstruct, so '
                          f'it cannot be silently raised for you. Got {self.numberOfThreads.get()} '
                          'threads.')

        # Reports the first missing tool only.
        for check in (Plugin.getProtacModelEnviron, Plugin.getProtacModelPython):
            try:
                check()
            except FileNotFoundError as e:
                errors.append(str(e))

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
    @staticmethod
    def _cleanReceptorOrTarget(inFile, outFile, ligandNameParam):
        """ Water-stripped copy of inFile keeping only the named warhead, if any, among
        heteroatoms. Upper-cased because PDB residue names are. """
        ligandName = ligandNameParam.get()
        if ligandName:
            cleanPDB(inFile, outFile, waters=True, hetatm=True,
                    het2keep=[ligandName.strip().upper()])
        else:
            cleanPDB(inFile, outFile, waters=True, hetatm=False)

    def _getSiteCoords(self):
        """ Parses siteCoords ("X,Y,Z") into a tuple of 3 floats, or None if malformed. """
        parts = self.siteCoords.get().strip().split(',')
        if len(parts) != 3:
            return None
        try:
            return tuple(float(p.strip()) for p in parts)
        except ValueError:
            return None

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
        if self.e3Ligand1.get() is not None and self.e3Ligand2.get() is not None:
            return 2
        return 1

    @staticmethod
    def _getSmiArg(smilesParam):
        """ 'none' is how PROTAC-Model means "no SMILES given". """
        value = smilesParam.get()
        return value.strip() if value else 'none'


