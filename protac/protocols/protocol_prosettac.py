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
from protac.constants import PROSETTAC_DIC, PROSETTAC_PYTHON_DIC, PROSETTAC_PYTHON2_DIC


class ProtPRosettaC(EMProtocol):
    """
    Models a PROTAC-mediated ternary complex with PRosettaC: relax, linker distance
    sampling, PatchDock global docking under that distance constraint, RosettaScripts
    local docking, constrained linker conformations and clustering of the best models.
    """
    _label = 'PROTAC ternary complex modeling (PRosettaC)'
    _devStatus = BETA
    _possibleOutputs = {'outputTernaryModels': SetOfAtomStructsChem}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        group = form.addGroup('Structures')
        group.addParam('structure1', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Structure 1 (static / PatchDock receptor)',
                       help='Static body for PatchDock. PRosettaC is ambiguous about '
                            'whether this should be the E3 ligase or the target, so the '
                            'fields are named by role. Do not pre-clean it: PRosettaC '
                            'extracts the chains itself, and renumbering the atoms would '
                            'break the anchor indices.')
        group.addParam('chain1', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 1 chain ID(s)',
                       help='One or more chain IDs as they appear in the PDB, e.g. "A" '
                            'or "AC".')
        group.addParam('structure2', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Structure 2 (moving / PatchDock ligand)',
                       help='Moving body for PatchDock. Do not pre-clean it either.')
        group.addParam('chain2', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 2 chain ID',
                       help='A single chain, different from those of Structure 1. The '
                            'final clustering is done on this chain.')

        group = form.addGroup('Warheads (Heads)')
        group.addParam('heads1', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       allowsNull=False, label='Structure 1 warhead set',
                       help='Set containing the warhead bound to Structure 1, in its '
                            'bound pose.')
        group.addParam('head1Name', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 1 warhead name',
                       help='Name of the molecule to use within the set above.')
        group.addParam('anchor1', params.IntParam, allowsNull=False,
                       label='Structure 1 anchor atom (1-based)',
                       help='1-based index, in the warhead file, of the atom where the '
                            'linker attaches. It must be a heavy atom. A wrong index gives '
                            'no clear error: the distance sampling just finds no '
                            'conformation.')
        group.addParam('heads2', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       allowsNull=False, label='Structure 2 warhead set',
                       help='Same as "Structure 1 warhead set", for Structure 2.')
        group.addParam('head2Name', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 2 warhead name', help='Same as "Structure 1 '
                            'warhead name", for this set.')
        group.addParam('anchor2', params.IntParam, allowsNull=False,
                       label='Structure 2 anchor atom (1-based)',
                       help='Same as "Structure 1 anchor atom", for this warhead.')

        group = form.addGroup('PROTAC')
        group.addParam('protacSmiles', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='PROTAC SMILES',
                       help='SMILES of the full PROTAC. Each warhead must be an exact '
                            'substructure of it, or no conformation is generated.')
        group.addParam('patchdockResults', params.IntParam, default=500,
                       label='PatchDock solutions to refine',
                       help='Best PatchDock solutions taken to local docking. The '
                            'published protocol uses 1000.')
        group.addParam('localNstruct', params.IntParam, default=10,
                       label='Local docking models per solution',
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
    def _insertAllSteps(self):
        # Scipion compares each step's arguments on its own when continuing a run, so every
        # step also gets a key with its predecessors' arguments: changing an input then
        # reruns the step that uses it and everything after it.
        chain1, chain2 = self.chain1.get().strip(), self.chain2.get().strip()
        stages = [
            (self.convertInputStep,
             [self.head1Name.get().strip(), self.head2Name.get().strip(),
              self.heads1.get().getObjId(), self.heads2.get().getObjId()]),
            (self.prepareStructuresStep,
             [self.structure1.get().getFileName(), chain1,
              self.structure2.get().getFileName(), chain2,
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

    def convertInputStep(self, head1Name, head2Name, heads1Id, heads2Id):
        """ Writes both warheads as SDF under the work dir. heads1Id/heads2Id are unused:
        they make Scipion rerun the step when the input sets change. """
        cleanPath(self._getWorkDir())
        os.makedirs(self._getWorkDir(), exist_ok=True)
        mol1 = self._findMolByName(self.heads1.get(), head1Name)
        mol2 = self._findMolByName(self.heads2.get(), head2Name)
        for mol, outFile in ((mol1, self._getHeadFile(1)), (mol2, self._getHeadFile(2))):
            srcFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
            convertToSdf(self, srcFile, sdfFile=outFile, overWrite=True)

    def prepareStructuresStep(self, struct1File, chain1, struct2File, chain2,
                              anchor1, anchor2, upstream):
        """ Phases 1+2 (entangled in the original, can't be split): per-structure
        addH_sdf -> translate_anchors -> clean -> mol_to_params -> relax -> clean. """
        Plugin.prepareRosettaScriptsShim(self._getRosettaShimDir())
        Plugin.preparePRosettaCBabelShim(self._getBabelShimDir())
        # Rosetta's clean_pdb.py only reads PDB; toPdb converts mmCIF without renumbering.
        struct1File = toPdb(struct1File, self._getExtraPath('structure1.pdb'))
        struct2File = toPdb(struct2File, self._getExtraPath('structure2.pdb'))

        # Absolute paths: the driver runs from the work dir.
        args = (f'--phase prepare '
               f'--struct1 "{os.path.abspath(struct1File)}" --chain1 "{chain1}" '
               f'--struct2 "{os.path.abspath(struct2File)}" --chain2 "{chain2}" '
               f'--head1 "{os.path.abspath(self._getHeadFile(1))}" --anchor1 {anchor1} '
               f'--head2 "{os.path.abspath(self._getHeadFile(2))}" --anchor2 {anchor2}')
        self._runDriver(args)

    def sampleDistStep(self, protacSmiles, upstream):
        """ Phase 3: linker distance sampling, PRosettaC's own pl.SampleDist(). """
        self._runDriver(f'--phase sampledist --smiles "{protacSmiles}"')

    def patchdockStep(self, globalResults, threshold, upstream):
        """ Phase 4: PatchDock global docking under the sampled distance constraint.
        Previous outputs are removed first: the results directory cannot exist yet. """
        cleanPath(self._getWorkDirFile('Patchdock_Results'))
        cleanPath(self._getWorkDirFile('Patchdock_params.txt'))

        self._runDriver(f'--phase patchdock --global-results {globalResults} '
                        f'--threshold {threshold}')

    def localDockingStep(self, nstruct, chain1, chain2, upstream):
        """ Phase 5: RosettaScripts local docking of every PatchDock solution. """
        self._runDriver(f'--phase localdocking --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--nstruct {nstruct} --threads {self.numberOfThreads.get()}')

    def constraintConfStep(self, chain1, chain2, upstream):
        """ Phase 6: constrained PROTAC conformations for each local docking solution. """
        self._runDriver(f'--phase constraintconf --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--threads {self.numberOfThreads.get()}')

    def clusteringStep(self, chain2, topScore, topLocal, rmsd, upstream):
        """ Phase 7: clustering of the best models by moving-chain RMSD. Previous
        outputs are removed first: clustering.main() fails if they exist. """
        cleanPath(self._getResultsDir())
        cleanPath(self._getWorkDirFile('result_summary.txt'))

        self._runDriver(f'--phase clustering --chain2 "{chain2}" --top-score {topScore} '
                        f'--top-local {topLocal} --cluster-rmsd {rmsd}')

    def createOutputStep(self, upstream):
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
        self._defineSourceRelation(self.structure1, outputSet)
        self._defineSourceRelation(self.structure2, outputSet)
        self._defineSourceRelation(self.heads1, outputSet)
        self._defineSourceRelation(self.heads2, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        chain1 = (self.chain1.get() or '').strip()
        chain2 = (self.chain2.get() or '').strip()
        if set(chain1) & set(chain2):
            errors.append('"Structure 1/2 chain ID(s)" must not overlap.')
        chains = chain1 + chain2
        if chains and (not chains.isalnum() or chains != chains.upper()
                       or set(chains) & set('XY')):
            errors.append('Chain IDs must be uppercase letters or digits, and cannot be X '
                          'or Y, which PRosettaC uses for the warheads.')
        # clustering.py compares line[21] == chain, so a multi-chain value selects no atom.
        if len(chain2) != 1:
            errors.append(f'"Structure 2 chain ID" must be a single chain. Got: "{chain2}".')

        if any((value or 0) < 1 for value in (self.patchdockResults.get(),
                                                self.localNstruct.get(),
                                                self.anchor1.get(), self.anchor2.get())):
            errors.append('Sampling sizes and anchor atoms (1-based) must be at least 1.')

        topScore, topLocal = self.clusterTopScore.get(), self.clusterTopLocal.get()
        if None not in (topScore, topLocal) and topLocal > topScore:
            errors.append('"Models to cluster by interface score" cannot be larger than '
                          '"Models to keep by total score".')

        for i, heads, name in ((1, self.heads1, self.head1Name), (2, self.heads2, self.head2Name)):
            if heads.get() is not None and name.get():
                try:
                    self._findMolByName(heads.get(), name.get().strip())
                except ValueError:
                    errors.append(f'"Structure {i} warhead name": no molecule named '
                                  f'"{name.get().strip()}" in that set.')
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
        Plugin.runCondaScript(Plugin.getPluginScript('run_prosettac.py'), args,
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

    @staticmethod
    def _findMolByName(smallMolSet, name):
        for mol in smallMolSet:
            if mol.getMolName() == name:
                return mol.clone()
        raise ValueError(f'No molecule named "{name}" found in {smallMolSet}.')

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
