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
Wrapper around PRosettaC (LondonLab/Karanicolas, https://github.com/LondonLab/PRosettaC),
an integrative modeling pipeline for PROTAC-mediated ternary complexes via PatchDock +
RosettaScripts.

The pipeline logic is not reimplemented: each step launches run_prosettac.py, which
calls PRosettaC's own functions. The only part replaced is its cluster layer: the jobs
the original submits to PBS/SGE/SLURM run in a local thread pool instead.
"""

import glob
import json
import os
import shutil
import csv


from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message, cleanPath
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pwem.convert.atom_struct import toPdb

from pwchem.objects import SetOfAtomStructsChem
from pwchem.utils import convertToSdf

from protac import Plugin
from protac.utils.molecules import findMolByName, listMolNames
from protac.constants import PROSETTAC_DIC, PROSETTAC_PYTHON_DIC, PROSETTAC_PYTHON2_DIC


class ProtPRosettaC(EMProtocol):
    """
    Models a PROTAC-mediated ternary complex with PRosettaC: relax, linker distance
    sampling, PatchDock global docking under that distance constraint, RosettaScripts
    local docking, constrained linker conformations and clustering of the best models.
    """
    _label = 'PRosettaC complex modeling'
    _devStatus = BETA
    _possibleOutputs = {'outputTernaryModels': SetOfAtomStructsChem}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        group = form.addGroup('Structures')
        group.addParam('e3ligase', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=True, label='E3 ligase structure: ',
                       help='Static body for PatchDock. PRosettaC is ambiguous about '
                            'whether this should be the E3 ligase or the target, so the '
                            'fields are named by role. Do not pre-clean it: PRosettaC '
                            'extracts the chains itself, and renumbering the atoms would '
                            'break the anchor indices.')
        group.addParam('targetProt', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=True, label='Target structure: ',
                       help='Moving body for PatchDock. Do not pre-clean it either.')

        group = form.addGroup('Warheads (Heads)')
        group.addParam('setE3Ligase', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       allowsNull=False, label='Set with E3 ligase ligand: ',
                       help='Set containing the warhead bound to Structure 1, in its '
                            'bound pose.')
        group.addParam('e3ligaseLigand', params.StringParam, allowsNull=False,
                       label='E3 ligase ligand: ',
                       help='Name of the molecule to use within the set above.')
        # todo get this with viewer
        group.addParam('anchor1', params.IntParam, allowsNull=False,
                       label='Anchor atom e3 ligase: ',
                       help='1-based index, in the warhead file, of the atom where the '
                            'linker attaches. It must be a heavy atom. A wrong index gives '
                            'no clear error: the distance sampling just finds no '
                            'conformation.')
        
        group.addParam('setTargetLigand', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       allowsNull=False, label='Set with target protein ligand: ',
                       help='Same as "Structure 1 warhead set", for Structure 2.')
        group.addParam('targetLIgand', params.StringParam, allowsNull=False,
                       label='Target ligand: ', help='Same as "Structure 1 '
                            'warhead name", for this set.')
        # todo get this with viewer
        group.addParam('anchor2', params.IntParam, allowsNull=False,
                       label='Anchor atom target: ',
                       help='Same as "Structure 1 anchor atom", for this warhead.')

        group = form.addGroup('PROTAC')
        group.addParam('protacMolecule', params.PointerParam,
                       pointerClass='SetOfSmallMolecules',
                       label='PROTAC molecule: ',
                       help='Full PROTAC molecule (E3 ligand - linker - '
                            'POI warhead), used to filter FRODOCK poses by compatibility '
                            'with the linker geometry.')

        group.addParam('patchdockResults', params.IntParam, default=500,
                       label='PatchDock solutions to refine: ',
                       help='Best PatchDock solutions taken to local docking. The '
                            'published protocol uses 1000.')
        group.addParam('localNstruct', params.IntParam, default=10,
                       label='Local docking models per solution: ',
                       help='Rosetta local docking models per PatchDock solution. The '
                            'published protocol uses 50.')

        group = form.addGroup('Advanced', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('patchdockThreshold', params.FloatParam, default=2.0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='PatchDock final clustering RMSD (A)',
                       help='Last field of PatchDock\'s "clusterParams" line. PRosettaC '
                            'uses 2.0 instead of the PatchDock default of 4.0. Lower values '
                            'keep more, more similar solutions.')
        group.addParam('clusterTopScore', params.IntParam, default=1000,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Models to keep by total score',
                       help='How many models to keep, ranked by total Rosetta score, '
                            'before clustering.')
        group.addParam('clusterTopLocal', params.IntParam, default=200,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Models to cluster by interface score',
                       help='How many of those models to cluster, ranked by interface '
                            'score.')
        group.addParam('clusterRmsd', params.FloatParam, default=4.0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Clustering RMSD threshold (A)',
                       help='DBSCAN radius, as RMSD over the CA atoms of the moving '
                            'chain.')

        form.addParallelSection(threads=4, mpi=1)

    @classmethod
    def _requireEnv(cls):
        for dic in (PROSETTAC_DIC, PROSETTAC_PYTHON_DIC, PROSETTAC_PYTHON2_DIC):
            Plugin._requireToolHome(dic)
        Plugin._requirePatchdockHome()
        Plugin.requireRosettaScriptsBinary()
        Plugin.requireObabelBinary()

    @classmethod
    def validateInstallation(cls):
        """ Replaces the plugin-wide check, which asks for PROTAC-Model's tools. """
        try:
            cls._requireEnv()
        except FileNotFoundError as e:
            return [str(e)]
        return []

    # --------------------------- STEPS functions ------------------------------
    def _insertAllStepsOld(self):
        # Scipion compares each step's arguments on its own when continuing a run, so every
        # step also gets a key with its predecessors' arguments: changing an input then
        # reruns the step that uses it and everything after it.
        chain1, chain2 = self.chain1.get().strip(), self.chain2.get().strip()
        stages = [
            (self.convertInputStep,
             [self.e3ligaseLigand.get().strip(), self.targetLIgand.get().strip(),
              self.setE3Ligase.get().getObjId(), self.setTargetLigand.get().getObjId()]),
            (self.prepareStructuresStep,
             [self.e3ligase.get().getFileName(), chain1,
              self.targetProt.get().getFileName(), chain2,
              self.anchor1.get(), self.anchor2.get()]),
            (self.sampleDistStep, [self.protacSmiles.get().strip()]),
            (self.patchdockStep, [self.patchdockResults.get(), self.patchdockThreshold.get()]),
            (self.localDockingStep, [self.localNstruct.get(), chain1, chain2]),
            (self.constraintConfStep, [chain1, chain2]),
            (self.clusteringStep, [chain2, self.clusterTopScore.get(),
                                   self.clusterTopLocal.get(), self.clusterRmsd.get()]),
            (self.createOutputStep, []),
        ]
        upstream, prerequisites = None, []
        for step, args in stages:
            key = [] if upstream is None else [json.dumps(upstream)]
            stepId = self._insertFunctionStep(step, *args, *key, prerequisites=prerequisites)
            upstream = (upstream or []) + args
            prerequisites = [stepId]

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.prepareStructuresStep)
        self._insertFunctionStep(self.sampleDistStep)
        self._insertFunctionStep(self.patchdockStep)
        self._insertFunctionStep(self.localDockingStep)
        self._insertFunctionStep(self.constraintConfStep)
        self._insertFunctionStep(self.clusteringStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertInputStep(self):
        """Copy the input structures to extra/, detect their chain IDs and
        rename conflicting target chains when necessary."""

        extraDir = self._getExtraPath()
        os.makedirs(extraDir, exist_ok=True)

        # 1. Copy/convert both structures to PDB files in extra/.
        # toPdb also handles mmCIF input.
        e3ligSrc = self.e3ligase.get().getFileName()
        targetSrc = self.targetProt.get().getFileName()

        e3ligFile = toPdb(
            e3ligSrc,
            self._getExtraPath('e3ligase.pdb')
        )
        targetFile = toPdb(
            targetSrc,
            self._getExtraPath('target.pdb')
        )

        # 2. Obtain all chain IDs from ATOM/HETATM records.
        e3Chains = self.getChainIds(e3ligFile)
        targetChains = self.getChainIds(targetFile)

        print(f'E3 ligase chains: {e3Chains}')
        print(f'Target chains: {targetChains}')

        # 3. Store the chain selections as strings for later steps.
        self.e3ChainIds = ''.join(e3Chains)
        self.targetChainIds = ''.join(targetChains)

        # 4. Rename only target chains that conflict with E3 chain IDs.
        conflicts = set(e3Chains) & set(targetChains)

        if conflicts:
            print(f'Conflicting chain IDs: {sorted(conflicts)}')

            # PDB chain IDs are single characters.
            # Exclude X/Y because PRosettaC reserves them for the warheads.
            availableIds = [
                c for c in 'ABCDEFGHIJKLMNOPQRSTUVWZ0123456789'
                if c not in e3Chains and c not in targetChains
            ]

            renameMap = {}

            for oldChain in targetChains:
                if oldChain in conflicts:
                    if not availableIds:
                        raise RuntimeError(
                            'Cannot rename conflicting target chains: '
                            'no unused PDB chain IDs remain.'
                        )

                    newChain = availableIds.pop(0)
                    renameMap[oldChain] = newChain

            # Rewrite the target PDB, modifying only the chain-ID column.
            tempFile = targetFile + '.tmp'

            with open(targetFile) as fin, open(tempFile, 'w') as fout:
                for line in fin:
                    if line.startswith(('ATOM  ', 'HETATM')):
                        oldChain = line[21].strip()

                        if oldChain in renameMap:
                            line = line[:21] + renameMap[oldChain] + line[22:]

                    fout.write(line)

            os.replace(tempFile, targetFile)

            print(f'Renamed target chains: {renameMap}')

            # Update the target chain string after renaming.
            self.targetChainIds = ''.join(
                renameMap.get(chain, chain) for chain in targetChains
            )

        # 5. Persist the chain selections for subsequent Scipion steps.
        chainInfo = {
            'e3ligase_file': os.path.abspath(e3ligFile),
            'target_file': os.path.abspath(targetFile),
            'e3_chains': self.e3ChainIds,
            'target_chains': self.targetChainIds,
        }

        with open(self._getExtraPath('chain_info.json'), 'w') as f:
            json.dump(chainInfo, f, indent=4)

        print(f'Final E3 ligase chains: {self.e3ChainIds}')
        print(f'Final target chains: {self.targetChainIds}')

        # 6. Obtain the SMILES of the full PROTAC molecule.
        protacSet = self.protacMolecule.get()
        mol = protacSet.getFirstItem()

        # Use the pose file for docked sets; otherwise, use the original file.
        if protacSet.isDocked():
            molFile = mol.getPoseFile()
        else:
            molFile = mol.getFileName()

        if not molFile:
            raise ValueError(
                f'No {"pose" if protacSet.isDocked() else "molecule"} '
                'file found for the PROTAC.'
            )

        molFile = os.path.abspath(molFile)
        ext = os.path.splitext(molFile)[1].lower()

        # Convert MOL files to SDF before RDKit processing.
        if ext == '.mol':
            sdfFile = os.path.join(
                extraDir, 'protac_ligand.sdf'
            )

            pwchemPlugin.runScript(
                self,
                'obabel_IO.py',
                [
                    '-i', molFile,
                    '-of', 'sdf',
                    '-o', 'protac_ligand',
                    '-od', os.path.join(extraDir,'prosettac')
                ],
                env=OPENBABEL_DIC,
                cwd=self._getPath()
            )

            if not os.path.isfile(sdfFile):
                raise FileNotFoundError(
                    f'SDF conversion failed for PROTAC: {sdfFile}'
                )

            molFile = sdfFile

        elif ext != '.sdf':
            raise ValueError(
                f'Unsupported PROTAC ligand format: {ext} ({molFile})'
            )

        # Generate the SMILES CSV using RDKit.
        smilesList = os.path.join(extraDir, 'ligands.txt')

        with open(smilesList, 'w') as f:
            f.write(f'{molFile}\n')

        smilesCsv = os.path.abspath(os.path.join(extraDir, 'ligand_smiles.csv'))

        pwchemPlugin.runScript(
            self,
            'rdkit_IO.py',
            [
                '-i', os.path.abspath(smilesList),
                '-of', 'smiles_csv',
                '-o', 'ligand_smiles',
                '-od', os.path.abspath(extraDir)
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

        if not rows or not rows[0].get(smilesColumn, '').strip():
            raise ValueError(
                f'No valid PROTAC SMILES found in {smilesCsv}.'
            )

        self.protacSmiles = rows[0][smilesColumn].strip()

        print(f'PROTAC SMILES: {self.protacSmiles}')

        # 7. Convert molecules to SDF
        workDir = self._getWorkDir()
        os.makedirs(workDir, exist_ok=True)

        ligands = [
            (self.getSpecifiedMol(
                self.setE3Ligase.get(),
                self.e3ligaseLigand.get()
            ), 'head1.sdf', self.setE3Ligase.get()),

            (self.getSpecifiedMol(
                self.setTargetLigand.get(),
                self.targetLIgand.get()
            ), 'head2.sdf', self.setTargetLigand.get()),

            (self.protacMolecule.get().getFirstItem(),
             'protac.sdf', self.protacMolecule.get())
        ]

        for mol, outName, molSet in ligands:
            # Use the pose file for docked sets; otherwise, use the original file.
            if molSet.isDocked():
                srcFile = mol.getPoseFile()
            else:
                srcFile = mol.getFileName()

            if not srcFile:
                raise ValueError(
                    f'No {"pose" if molSet.isDocked() else "molecule"} '
                    f'file found for {outName}.'
                )

            outFile = os.path.abspath(
                os.path.join(workDir, outName)
            )

            convertToSdf(
                self,
                srcFile,
                sdfFile=outFile,
                overWrite=True
            )

            if not os.path.isfile(outFile) or os.path.getsize(outFile) == 0:
                raise FileNotFoundError(
                    f'SDF conversion failed for {outName}: {outFile}'
                )

            print(f'Converted {srcFile} -> {outFile}')


    def _getChainInfo(self):
        with open(self._getExtraPath('chain_info.json')) as f:
            return json.load(f)


    def convertInputStepOld(self, e3ligaseLigand, targetLIgand, setE3LigaseId, setTargetLigandId):
        """ Writes both warheads as SDF under the work dir. setE3LigaseId/setTargetLigandId are unused:
        they make Scipion rerun the step when the input sets change. """
        cleanPath(self._getWorkDir())
        os.makedirs(self._getWorkDir(), exist_ok=True)
        mol1 = self._findMolByName(self.setE3Ligase.get(), e3ligaseLigand)
        mol2 = self._findMolByName(self.setTargetLigand.get(), targetLIgand)
        for mol, outFile in ((mol1, self._getHeadFile(1)), (mol2, self._getHeadFile(2))):
            srcFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
            convertToSdf(self, srcFile, sdfFile=outFile, overWrite=True)

    def prepareStructuresStep(self):
        """ Phases 1+2 (entangled in the original, can't be split): per-structure
        addH_sdf -> translate_anchors -> clean -> mol_to_params -> relax -> clean. """
        chainInfo = self._getChainInfo()

        struct1File = chainInfo["e3ligase_file"]
        chain1 = chainInfo["e3_chains"]
        struct2File = chainInfo["target_file"]
        chain2 = chainInfo["target_chains"]
        anchor1 = self.anchor1.get()
        anchor2 = self.anchor2.get()

        Plugin.prepareRosettaScriptsShim(self._getRosettaShimDir())
        Plugin.preparePRosettaCBabelShim(self._getBabelShimDir())
        # Rosetta's clean_pdb.py only reads PDB; toPdb converts mmCIF without renumbering.
        struct1File = toPdb(struct1File, self._getExtraPath('e3ligase.pdb'))
        struct2File = toPdb(struct2File, self._getExtraPath('targetProt.pdb'))

        # Absolute paths: the driver runs from the work dir.
        args = (f'--phase prepare '
               f'--struct1 "{os.path.abspath(struct1File)}" --chain1 "{chain1}" '
               f'--struct2 "{os.path.abspath(struct2File)}" --chain2 "{chain2}" '
               f'--head1 "{os.path.abspath(self._getHeadFile(1))}" --anchor1 {anchor1} '
               f'--head2 "{os.path.abspath(self._getHeadFile(2))}" --anchor2 {anchor2}')
        self._runDriver(args)

    def sampleDistStep(self):
        """ Phase 3: linker distance sampling, PRosettaC's own pl.SampleDist(). """
        smilesCsv = self._getExtraPath('ligand_smiles.csv')
        if not os.path.isfile(smilesCsv):
            raise FileNotFoundError(
                f'PROTAC SMILES CSV not found: {smilesCsv}'
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
        if not rows or not rows[0].get(smilesColumn, '').strip():
            raise ValueError(
                f'No valid PROTAC SMILES found in {smilesCsv}.'
            )

        protacSmiles = rows[0][smilesColumn].strip()
        self._runDriver(f'--phase sampledist --smiles "{protacSmiles}"')

    def patchdockStep(self):
        """ Phase 4: PatchDock global docking under the sampled distance constraint.
        Previous outputs are removed first: the results directory cannot exist yet. """
        globalResults = self.patchdockResults.get()
        threshold = self.patchdockThreshold.get()

        cleanPath(self._getWorkDirFile('Patchdock_Results'))
        cleanPath(self._getWorkDirFile('Patchdock_params.txt'))

        self._runDriver(f'--phase patchdock --global-results {globalResults} '
                        f'--threshold {threshold}')

    def localDockingStep(self):
        """ Phase 5: RosettaScripts local docking of every PatchDock solution. """
        nstruct = self.localNstruct.get()
        chainInfo = self._getChainInfo()
        chain1 = chainInfo["e3_chains"]
        chain2 = chainInfo["target_chains"]
        self._runDriver(f'--phase localdocking --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--nstruct {nstruct} --threads {self.numberOfThreads.get()}')

    def constraintConfStep(self):
        """ Phase 6: constrained PROTAC conformations for each local docking solution. """
        chainInfo = self._getChainInfo()
        chain1 = chainInfo["e3_chains"]
        chain2 = chainInfo["target_chains"]

        self._runDriver(f'--phase constraintconf --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--threads {self.numberOfThreads.get()}')

    def clusteringStep(self):
        """ Phase 7: clustering of the best models by moving-chain RMSD. Previous
        outputs are removed first: clustering.main() fails if they exist. """
        chainInfo = self._getChainInfo()
        topScore = self.clusterTopScore.get()
        topLocal = self.clusterTopLocal.get()
        rmsd = self.clusterRmsd.get()
        chain2 = chainInfo["target_chains"]
        cleanPath(self._getResultsDir())
        cleanPath(self._getWorkDirFile('result_summary.txt'))

        self._runDriver(f'--phase clustering --chain2 "{chain2}" --top-score {topScore} '
                        f'--top-local {topLocal} --cluster-rmsd {rmsd}')

    def createOutputStep(self):
        """ One AtomStruct per clustered model, with its cluster (1 is the best ranked)
        and its total Rosetta score. """
        clusterDirs = self._getClusterDirs()
        if not clusterDirs:
            raise RuntimeError(
                f'{self._getResultsDir()} holds no cluster. PRosettaC clustered no model: '
                'either no local docking solution admitted a valid conformation of the '
                'full PROTAC (check the warhead SDFs, their anchor atoms and the PROTAC '
                'SMILES), or none of the resulting models scored below the 0 energy '
                'threshold.')

        scores = self._parseScoreFile(
            os.path.join(self._getWorkDirFile('Patchdock_Results'), 'score.sc'))

        outputSet = SetOfAtomStructsChem().create(self._getPath())
        for clusterId, clusterDir in clusterDirs:
            for modelFile in sorted(glob.glob(os.path.join(clusterDir, '*.pdb'))):
                modelName = os.path.splitext(os.path.basename(modelFile))[0]
                atomStruct = AtomStruct(filename=modelFile)
                atomStruct.setObjLabel(f'cluster{clusterId}_{modelName}')
                # Set on every item, even without a score, so all items share attributes.
                atomStruct._clusterId = pwobj.Integer(clusterId)
                atomStruct._score = pwobj.Float(scores.get(modelName))
                outputSet.append(atomStruct)

        if len(outputSet) == 0:
            raise RuntimeError(
                f'{self._getResultsDir()} has cluster directories but no model .pdb in '
                'them. This should not happen - check the clustering step log.')

        self._defineOutputs(outputTernaryModels=outputSet)
        self._defineSourceRelation(self.e3ligase, outputSet)
        self._defineSourceRelation(self.targetProt, outputSet)
        self._defineSourceRelation(self.setE3Ligase, outputSet)
        self._defineSourceRelation(self.setTargetLigand, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []

        return errors

    def _summary(self):
        summary = []
        if self.isFinished() and hasattr(self, 'outputTernaryModels'):
            models = self.outputTernaryModels
            if len(models):
                clusters = {model._clusterId.get() for model in models}
                summary.append(f'{len(models)} ternary complex model(s) in '
                               f'{len(clusters)} cluster(s). Cluster 1 is the best '
                               'ranked: clusters are ordered by size first, then by '
                               'average Rosetta score.')
            summaryFile = self._getWorkDirFile('result_summary.txt')
            if os.path.exists(summaryFile):
                with open(summaryFile) as f:
                    summary.extend(line.strip() for line in f if line.strip())
        return summary

    def _citations(self):
        return ['Zaidman2020']

    # --------------------------- UTILS functions ------------------------------
    def _runDriver(self, args):
        """ Runs one phase of run_prosettac.py from the work dir. """
        Plugin.runCondaScript(self, Plugin.getPluginScript('run_prosettac.py'), args,
                              PROSETTAC_PYTHON_DIC,
                              extraEnvDict=Plugin.getPRosettaCEnviron(
                                  self._getRosettaShimDir(), self._getBabelShimDir()),
                              cwd=self._getWorkDir())

    def _getClusterDirs(self):
        """ [(clusterId, path)] for every Results/cluster<N>, sorted numerically by N. """
        clusters = []
        for path in glob.glob(os.path.join(self._getResultsDir(), 'cluster*')):
            suffix = os.path.basename(path)[len('cluster'):]
            if os.path.isdir(path) and suffix.isdigit():
                clusters.append((int(suffix), path))
        return sorted(clusters)

    @staticmethod
    def _parseScoreFile(scoreFile):
        """ Rosetta score file -> {description: total_score}, looking columns up by name.
        Returns {} if the file is missing. """
        scores = {}
        if not os.path.exists(scoreFile):
            return scores
        columns = None
        with open(scoreFile) as f:
            for line in f:
                fields = line.split()
                if not fields or fields[0] != 'SCORE:':
                    continue
                if columns is None:
                    columns = {name: i for i, name in enumerate(fields)}
                    continue
                try:
                    scores[fields[columns['description']]] = \
                        float(fields[columns['total_score']])
                except (KeyError, IndexError, ValueError):
                    continue
        return scores

    def gete3ligaseLigands(self):
        """ Molecule names in the Structure 1 warhead set, listed by its name wizard. """
        return listMolNames(self.setE3Ligase.get())

    def gettargetLIgands(self):
        """ See gete3ligaseLigands. """
        return listMolNames(self.setTargetLigand.get())

    @staticmethod
    def _findMolByName(smallMolSet, name):
        return findMolByName(smallMolSet, name)

    def _getWorkDir(self):
        return self._getExtraPath('prosettac')

    def _getWorkDirFile(self, name):
        return self._getExtraPath('prosettac', name)

    def _getHeadFile(self, i):
        return self._getExtraPath('prosettac', f'head{i}.sdf')

    def _getResultsDir(self):
        return self._getExtraPath('prosettac', 'Results')

    def _getRosettaShimDir(self):
        return self._getExtraPath('rosetta_shim')

    def _getBabelShimDir(self):
        return self._getExtraPath('babel_shim')

    def getSpecifiedMol(self, inputSmallMolecules, inputLigand):
        myMol = None

        for mol in inputSmallMolecules:
            if mol.__str__() == inputLigand:
                myMol = mol.clone()
                break

        if myMol is None:
            raise ValueError(
                f'The input ligand "{inputLigand}" was not found '
                'in the selected molecule set.'
            )

        return myMol

    def getChainIds(self, pdbFile):
        chains = set()

        with open(pdbFile) as f:
            for line in f:
                if line.startswith(('ATOM  ', 'HETATM')):
                    chain = line[21].strip()
                    if chain:
                        chains.add(chain)

        return sorted(chains)
