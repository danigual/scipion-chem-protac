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

All 7 phases of the original pipeline are covered: structure/warhead preparation and
relax, linker distance sampling, PatchDock global docking, RosettaScripts local docking,
constrained linker conformations, and clustering of the top models.

None of PRosettaC's own pipeline logic is reimplemented here: each step stages arguments
and launches protac/scripts/run_prosettac.py, which calls straight into PRosettaC's own
utils/rosetta/protac_lib/clustering functions. The one thing the driver replaces rather
than calls is PRosettaC's cluster layer: the original submits phases 5 and 6 to a
PBS/SGE/SLURM scheduler, and the driver runs the very same commands in a local thread
pool instead (see run_prosettac.py's _runCommandsInParallel).
"""

import glob
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
    Models a PROTAC-mediated ternary complex (two proteins bridged by a PROTAC) via
    PRosettaC: structure/warhead preparation and relax, sampling of the distance between
    the two PROTAC anchor points, PatchDock global docking under that distance
    constraint, RosettaScripts local docking of each PatchDock solution, generation of
    constrained conformations of the full linker for each of them, and clustering of the
    top-scoring models into ranked ternary complex predictions.
    """
    _label = 'PROTAC ternary complex modeling (PRosettaC)'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        group = form.addGroup('Structures')
        group.addParam('structure1', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Structure 1 (static / PatchDock receptor)',
                       help="PRosettaC's own README and its clustering.py contradict each "
                            'other on whether this should be the E3 ligase or the target '
                            'protein - the only thing the code actually requires is that '
                            "this is the PatchDock receptor (static). Don't clean it with "
                            "a separate cleanPDB step first: PRosettaC's own rs.clean() "
                            'does the chain extraction, and pre-cleaning would renumber '
                            'atoms and break the anchor atom indices below.')
        group.addParam('chain1', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 1 chain ID(s)',
                       help='One or more chain IDs, as they appear in the PDB (e.g. "A" '
                            'or "AC"). Passed as-is to PRosettaC.')
        group.addParam('structure2', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Structure 2 (moving / PatchDock ligand)',
                       help='The other protein - PatchDock ligand (moving body). Same '
                            'caveat about not pre-cleaning as Structure 1.')
        group.addParam('chain2', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 2 chain ID', help='Same as "Structure 1 chain '
                            'ID(s)", for this structure. Must not share any letter with '
                            'it, and must be a single chain: the final clustering is done '
                            'on this structure and PRosettaC matches the PDB chain column '
                            'against this value one character at a time.')

        group = form.addGroup('Warheads (Heads)')
        group.addParam('heads1', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       allowsNull=False, label='Structure 1 warhead set',
                       help='Set containing the warhead already bound to Structure 1, in '
                            'its bound conformation/orientation.')
        group.addParam('head1Name', params.StringParam, allowsNull=False,
                       validators=[params.NonEmptyCondition()],
                       label='Structure 1 warhead name',
                       help='Name of the molecule to use within the set above.')
        group.addParam('anchor1', params.IntParam, allowsNull=False,
                       label='Structure 1 anchor atom (1-based)',
                       help="1-based atom index, within the warhead file as given, of "
                            "PRosettaC's anchor atom. Must be a heavy atom, not a hydrogen "
                            '(hydrogens are stripped before use and the index is remapped). '
                            'The warhead must be 3D. Manual only: a wrong index does not '
                            'raise a clear error, it surfaces later as SampleDist finding '
                            'no PROTAC conformation to sample (0,0).')
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
                       label='PROTAC SMILES', help='SMILES of the full PROTAC molecule.')
        group.addParam('doFull', params.BooleanParam, default=False,
                       label='Full run', help='Sampling depth, as in PRosettaC itself: '
                            '1000 vs 500 PatchDock solutions to take forward, and 50 vs '
                            '10 Rosetta local docking models per solution. A full run is '
                            'the published protocol; a quick one is roughly 10x cheaper.')

        group = form.addGroup('Advanced', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('patchdockThreshold', params.FloatParam, default=2.0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='PatchDock final clustering RMSD (A)',
                       help='Last field of PatchDock\'s own "clusterParams" line, which '
                            'PRosettaC overrides from PatchDock\'s 4.0 default down to '
                            '2.0. Lowering it keeps more, more similar solutions.')
        group.addParam('clusterTopScore', params.IntParam, default=1000,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Models to keep by total score',
                       help='First stage of the final clustering: how many models to keep, '
                            'ranked by their total Rosetta score. 1000 is the published '
                            'value.')
        group.addParam('clusterTopLocal', params.IntParam, default=200,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Models to cluster by interface score',
                       help='Second stage: how many of the models kept above to actually '
                            'cluster, ranked this time by interface score. 200 is the '
                            'published value.')
        group.addParam('clusterRmsd', params.FloatParam, default=4.0,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Clustering RMSD threshold (A)',
                       help='DBSCAN neighbourhood radius: two models closer than this '
                            '(RMSD over the moving chain CA atoms) fall in the same '
                            'cluster. 4.0 A is the published value.')

        form.addParallelSection(threads=4, mpi=1)

    @classmethod
    def _requireEnv(cls):
        for dic in (PROSETTAC_DIC, PROSETTAC_PYTHON_DIC, PROSETTAC_PYTHON2_DIC):
            Plugin._requireToolHome(dic)
        Plugin._requirePatchdockHome()
        # Not just "the variable is set": this also checks the directory exists and holds
        # a rosetta_scripts binary, which is the part PRosettaC actually needs and the
        # only one that fails late (in the middle of prepareStructuresStep) otherwise.
        Plugin.requireRosettaScriptsBinary()
        Plugin.requireObabelBinary()

    @classmethod
    def validateInstallation(cls):
        """ Overrides the plugin-wide default (Plugin.validateInstallation, checked
        against _pathVars - PROTAC_MODEL_HOME, a dependency this protocol never touches),
        same pattern as ProtPROTACTransplantWarhead. A classmethod because that is what
        Protocol.isInstalled() calls it as. _validate() does not call it: Protocol.validate()
        already collects it separately, and calling it again duplicates every error. """
        try:
            cls._requireEnv()
        except FileNotFoundError as e:
            return [str(e)]
        return []

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        convertId = self._insertFunctionStep(
            self.convertInputStep, self.head1Name.get().strip(), self.head2Name.get().strip(),
            self.heads1.get().getObjId(), self.heads2.get().getObjId(),
            prerequisites=[])
        prepareId = self._insertFunctionStep(
            self.prepareStructuresStep,
            self.structure1.get().getFileName(), self.chain1.get().strip(),
            self.structure2.get().getFileName(), self.chain2.get().strip(),
            self.anchor1.get(), self.anchor2.get(),
            prerequisites=[convertId])
        sampleDistId = self._insertFunctionStep(
            self.sampleDistStep, self.protacSmiles.get().strip(), prerequisites=[prepareId])
        patchdockId = self._insertFunctionStep(
            self.patchdockStep, self.doFull.get(), self.patchdockThreshold.get(),
            prerequisites=[sampleDistId])
        localDockId = self._insertFunctionStep(
            self.localDockingStep, self.doFull.get(), self.chain1.get().strip(),
            self.chain2.get().strip(), prerequisites=[patchdockId])
        constraintId = self._insertFunctionStep(
            self.constraintConfStep, self.chain1.get().strip(), self.chain2.get().strip(),
            prerequisites=[localDockId])
        clusteringId = self._insertFunctionStep(
            self.clusteringStep, self.chain2.get().strip(), self.clusterTopScore.get(),
            self.clusterTopLocal.get(), self.clusterRmsd.get(),
            prerequisites=[constraintId])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[clusteringId])

    def convertInputStep(self, head1Name, head2Name, heads1Id, heads2Id):
        """ Resolves the two warhead SDFs into real .sdf files under the work dir
        (SampleDist/translate_anchors need Chem.SDMolSupplier). Doesn't touch structure1/
        structure2 - see their help text for why.

        heads1Id/heads2Id are not read here: they are in the signature so that the step's
        argument string, which is what Scipion compares to decide whether to re-execute a
        step, changes when the input set does. Without them, swapping in a different
        SetOfSmallMolecules holding a molecule of the same name would leave this step
        green with the old warhead. """
        cleanPath(self._getWorkDir())
        os.makedirs(self._getWorkDir(), exist_ok=True)
        mol1 = self._findMolByName(self.heads1.get(), head1Name)
        mol2 = self._findMolByName(self.heads2.get(), head2Name)
        for mol, outFile in ((mol1, self._getHeadFile(1)), (mol2, self._getHeadFile(2))):
            srcFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
            convertToSdf(self, srcFile, sdfFile=outFile, overWrite=True)

    def prepareStructuresStep(self, struct1File, chain1, struct2File, chain2,
                              anchor1, anchor2):
        """ Phases 1+2 (entangled in the original, can't be split): per-structure
        addH_sdf -> translate_anchors -> clean -> mol_to_params -> relax -> clean. """
        rosettaShim = Plugin.prepareRosettaScriptsShim(self._getRosettaShimDir())
        babelShim = Plugin.preparePRosettaCBabelShim(self._getBabelShimDir())
        # Rosetta's clean_pdb.py only reads PDB; toPdb converts mmCIF without renumbering.
        struct1File = toPdb(struct1File, self._getExtraPath('structure1.pdb'))
        struct2File = toPdb(struct2File, self._getExtraPath('structure2.pdb'))

        # Absolute + quoted: the driver runs with cwd=self._getWorkDir(), so a relative
        # path here would resolve against the wrong directory (same reasoning as
        # frodockStep in protocol_protac_model.py).
        args = (f'--phase prepare '
               f'--struct1 "{os.path.abspath(struct1File)}" --chain1 "{chain1}" '
               f'--struct2 "{os.path.abspath(struct2File)}" --chain2 "{chain2}" '
               f'--head1 "{os.path.abspath(self._getHeadFile(1))}" --anchor1 {anchor1} '
               f'--head2 "{os.path.abspath(self._getHeadFile(2))}" --anchor2 {anchor2}')
        self._runDriver(args, rosettaShim=rosettaShim, babelShim=babelShim)

    def sampleDistStep(self, protacSmiles):
        """ Phase 3: linker distance sampling, PRosettaC's own pl.SampleDist(). """
        self._runDriver(f'--phase sampledist --smiles "{protacSmiles}"')

    def patchdockStep(self, doFull, threshold):
        """ Phase 4: PatchDock global docking under the sampled distance constraint,
        PRosettaC's own utils.patchdock(). That function does 'os.mkdir' without exist_ok
        and writes fixed top-level filenames, so a rerun needs these cleared first or it
        crashes on FileExistsError instead of redoing the docking. """
        cleanPath(self._getWorkDirFile('Patchdock_Results'))
        cleanPath(self._getWorkDirFile('Patchdock_cst'))
        cleanPath(self._getWorkDirFile('Patchdock_params.txt'))

        globalResults = 1000 if doFull else 500
        self._runDriver(f'--phase patchdock --global-results {globalResults} '
                        f'--threshold {threshold}')

    def localDockingStep(self, doFull, chain1, chain2):
        """ Phase 5: RosettaScripts local docking of every PatchDock solution, one job per
        solution, each command line built by PRosettaC's own rs.local_docking(). The
        original submits those jobs to a PBS/SGE/SLURM scheduler; the driver runs the very
        same commands locally, at most "Threads" at a time. Files left over from a
        previous attempt are cleared by the driver itself, next to where they are written
        (Patchdock_Results/). """
        nstruct = 50 if doFull else 10
        self._runDriver(f'--phase localdocking --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--nstruct {nstruct} --threads {self.numberOfThreads.get()}')

    def constraintConfStep(self, chain1, chain2):
        """ Phase 6: constrained conformations of the full linker for each local docking
        solution, one PRosettaC constraint_generation.py job per solution. Same
        scheduler-vs-local-pool substitution as localDockingStep, and same reason for the
        cleanup living in the driver. """
        self._runDriver(f'--phase constraintconf --chain1 "{chain1}" --chain2 "{chain2}" '
                        f'--threads {self.numberOfThreads.get()}')

    def clusteringStep(self, chain2, topScore, topLocal, rmsd):
        """ Phase 7: DBSCAN clustering of the top models by moving-chain RMSD, PRosettaC's
        own clustering.main(). It does 'os.mkdir(Results/)' without exist_ok and writes
        result_summary.txt under a fixed name, so a rerun needs both cleared or it dies on
        FileExistsError instead of reclustering. """
        cleanPath(self._getResultsDir())
        cleanPath(self._getWorkDirFile('result_summary.txt'))

        self._runDriver(f'--phase clustering --chain2 "{chain2}" --top-score {topScore} '
                        f'--top-local {topLocal} --cluster-rmsd {rmsd}')

    def createOutputStep(self):
        """ Collects the clustered ternary complex models into a SetOfAtomStructsChem, one
        AtomStruct per model, tagged with its cluster (cluster 1 is the best ranked -
        clustering.py ranks clusters by size first, then by average score) and with its
        total Rosetta score.

        Two caveats about that score. It is the *total* score of the model (score.sc),
        which is what the cluster ranking averages - not the interface score (I_sc, in
        local.fasc) that clustering.py uses one step earlier to pick which models to
        cluster at all; recovering I_sc per model means inverting clustering.py's own
        model-name mangling, which is not done here. And despite the shared attribute
        name, it is not comparable to ProtPROTACModel's _score, which is a VoroMQA
        interface score: same name, different quantity and different scale. """
        clusterDirs = self._getClusterDirs()
        # Without this, a pipeline that modelled nothing would fall through to an empty
        # output set and a green protocol - worse than a clear failure, since nothing
        # downstream would flag it.
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
                # Dynamic attributes (AtomStruct has no score/cluster field), same
                # mechanism ProtPROTACModel uses for its own per-item score. Both are set
                # on every item, even when the score is unknown (None), so that every
                # element of the set carries the same attributes.
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
        overlap = set(self.chain1.get().strip()) & set(self.chain2.get().strip())
        if overlap:
            errors.append('"Structure 1/2 chain ID(s)" must not overlap (shared: '
                          f'{"".join(sorted(overlap))}).')

        # PRosettaC's clustering.py matches the PDB chain column one character at a time
        # ("line[21] == chain"), so a multi-chain moving body selects no CA atom at all
        # and it dies on a division by zero - after every docking phase has already run.
        if len(self.chain2.get().strip()) != 1:
            errors.append('"Structure 2 chain ID(s)" must be a single chain: PRosettaC '
                          'clusters the final models on the moving structure and compares '
                          'the PDB chain column against this value character by character. '
                          f'Got: "{self.chain2.get().strip()}".')

        return errors

    def _summary(self):
        summary = []
        # Set/Object.get() is for scalar attributes and always returns None on a Set -
        # isFinished() plus the hasattr check (set by _defineOutputs) is the right test,
        # same as ProtPROTACModel's own _summary().
        if self.isFinished() and hasattr(self, 'outputTernaryModels'):
            models = self.outputTernaryModels
            if len(models):
                clusters = {model._clusterId.get() for model in models}
                summary.append(f'{len(models)} ternary complex model(s) in '
                               f'{len(clusters)} cluster(s). Cluster 1 is the best '
                               'ranked: clusters are ordered by size first, then by '
                               'average Rosetta score.')
            # PRosettaC writes its own run-wide counts (local docking solutions, models
            # below the energy threshold, clusters with at least 5 members) - worth
            # surfacing verbatim rather than recomputing.
            summaryFile = self._getWorkDirFile('result_summary.txt')
            if os.path.exists(summaryFile):
                with open(summaryFile) as f:
                    summary.extend(line.strip() for line in f if line.strip())
        return summary

    def _citations(self):
        # protac/bibtex.py doesn't exist yet, any key here would be
        # silently dropped - see ProtPROTACTransplantWarhead's own _citations.
        return []

    # --------------------------- UTILS functions ------------------------------
    def _runDriver(self, args, rosettaShim=None, babelShim=None):
        """ Launches run_prosettac.py for one phase. Always from the work dir, and always
        with the full environment: PRosettaC's utils.py reads its 4 variables at import
        time, so every phase needs all of them regardless of which it uses.
        rosettaShim/babelShim: only prepareStructuresStep passes them, since it is the
        step that builds the shims; everywhere else the same paths are recomputed. """
        Plugin.runCondaScript(Plugin.getPluginScript('run_prosettac.py'), args,
                              PROSETTAC_PYTHON_DIC,
                              extraEnvDict=Plugin.getPRosettaCEnviron(
                                  rosettaHome=rosettaShim or self._getRosettaShimDir(),
                                  obDir=babelShim or self._getBabelShimDir()),
                              cwd=self._getWorkDir())

    def _getClusterDirs(self):
        """ [(clusterId, path)] for every Results/cluster<N>, ordered by N - which is the
        ranking clustering.py assigned them (1 = best). Sorted numerically on purpose:
        a plain sort would put cluster10 before cluster2. """
        clusters = []
        for path in glob.glob(os.path.join(self._getResultsDir(), 'cluster*')):
            suffix = os.path.basename(path)[len('cluster'):]
            if os.path.isdir(path) and suffix.isdigit():
                clusters.append((int(suffix), path))
        return sorted(clusters)

    @staticmethod
    def _parseScoreFile(scoreFile):
        """ A Rosetta score file -> {model description: total score}. Columns are looked
        up by name in the header line rather than by position (which is what PRosettaC's
        own clustering.py does) so that a Rosetta build with extra score terms cannot
        silently shift them. Returns {} if the file is missing or has no header: the
        score is metadata on the output, not a reason to fail a finished run. """
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
        """ Where clustering.py writes its ranked clusters (one cluster<N>/ per cluster). """
        return self._getExtraPath('prosettac', 'Results')

    def _getRosettaShimDir(self):
        return self._getExtraPath('rosetta_shim')

    def _getBabelShimDir(self):
        return self._getExtraPath('babel_shim')
