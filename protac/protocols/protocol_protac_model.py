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

None of PROTAC-Model's own pipeline logic is reimplemented here: each phase below just
stages arguments and launches protac/scripts/run_protac_model.py (a Python 2 driver,
run under a dedicated conda env), which in turn calls straight into PROTAC-Model's own
utils.frodock/utils.rosetta functions.
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
                            'This matters more than it looks: PROTAC-Model reads *every* '
                            'HETATM record as part of the warhead, without filtering by '
                            'residue. Leaving a sulfate or a glycerol from the '
                            'crystallization buffer in the file merges it into the ligand, '
                            'which shifts the anchoring point and yields wrong ternary models '
                            'with perfectly normal-looking scores.\n\n'
                            'Leave empty only if the structure already contains the warhead '
                            'as its single heteroatom.')
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
                       help='Coordinates of a point on the receptor surface, at the '
                            'interface where the ternary complex is expected to form. '
                            'FRODOCK uses this point to restrict the global docking search.\n\n'
                            "It must be the centroid of the receptor's warhead HETATM "
                            'records, computed on the same structure PROTAC-Model will '
                            'actually see: if "Receptor warhead residue name" is set '
                            'above, that means the *cleaned* structure (warhead only), '
                            'not the raw input - the two can differ once other '
                            'heteroatoms are stripped out.')

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
        """ Converts the receptor/target AtomStructs into the clean PDBs that FRODOCK
        expects. Water is always stripped. Other heteroatoms are stripped down to just
        the named warhead when its residue name is given (receptorLigandName/
        targetLigandName), since PROTAC-Model itself copies *every* HETATM record into
        the ligand file without filtering by residue - left unfiltered, a stray sulfate
        or glycerol from the crystallization buffer would merge into the warhead and
        shift its anchoring point. Left empty, the field is a no-op: every non-water
        heteroatom is kept, which is only safe if the structure already has the warhead
        as its sole heteroatom. """
        self._cleanReceptorOrTarget(self.inputReceptor.get().getFileName(),
                                    self._getReceptorFile(), self.receptorLigandName)
        self._cleanReceptorOrTarget(self.inputTarget.get().getFileName(),
                                    self._getTargetFile(), self.targetLigandName)

        # Built once here, reused by frodockStep/filterPosesStep/refineStep (see
        # _getFrodockBinDir): picks whichever FRODOCK build (intel/gcc) actually loads on
        # this machine, since PROTAC-Model's own fallback only checks the intel build's
        # file exists, not that it runs (see Plugin.prepareFrodockBinDir).
        Plugin.prepareFrodockBinDir(self._getFrodockBinDir())

    def frodockStep(self):
        """ Runs FRODOCK global docking by launching run_protac_model.py --phase frodock,
        which stages its inputs (chain renaming, protac.smi, optional E3 ligand SDFs) and
        then calls PROTAC-Model's own fro.frodock(site) - none of that logic is
        reimplemented here. """
        # extra/frodock/: the driver's --phase frodock/filter both run from this cwd
        # (PROTAC-Model's own fro.frodock()/fro.filter_frodock() read/write files
        # relative to it, e.g. 'receptor.pdb', 'protac.smi').
        frodockDir = self._getExtraPath('frodock')
        os.makedirs(frodockDir, exist_ok=True)

        # Absolute + quoted paths. Absolute because _getExtraPath()/_getPath() return
        # paths relative to the Scipion project directory, while every process launched
        # below runs with a cwd of our own choosing (frodockDir here, rosettaDir in
        # refineStep) - a relative path would resolve against the wrong directory inside
        # the driver. Quoted because the underlying command goes through a shell, which
        # would otherwise split a spaced path into extra tokens.
        # This applies to CLI arguments only: paths that travel as environment variables,
        # as the driver script path, or as cwd are made absolute centrally instead (see
        # Plugin.getProtacModelEnviron/runCondaScript/runProgram), so they cannot be
        # forgotten at a call site.
        receptorFile = os.path.abspath(self._getReceptorFile())
        targetFile = os.path.abspath(self._getTargetFile())
        # From the parsed tuple, not the raw form string, so stray whitespace can't split
        # --site into extra shell tokens.
        x, y, z = self._getSiteCoords()
        site = f'{x:.4f},{y:.4f},{z:.4f}'
        args = (f'--phase frodock --receptor "{receptorFile}" --target "{targetFile}" '
               f'--smiles "{self.protacSmiles.get().strip()}" --site={site}')

        # Both E3 ligand conformers, or neither: --phase filter derives lig_locate_num
        # from this same pair (see _getLigLocateNum) - _validate() already rejects a lone
        # conformer.
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
        """ Filters the FRODOCK poses by compatibility with the PROTAC/warhead geometry,
        via run_protac_model.py --phase filter -> PROTAC-Model's own fro.filter_frodock().
        Same cwd as frodockStep: filter_frodock() reads the files frodock() just wrote
        there (frodock_score.txt, receptor.pdb, target.pdb...). """
        # PROTAC-Model's filter_frodock() only creates extra/frodock_results if it is
        # missing (utils/frodock.py: 'if os.path.exists(...) == False: os.makedirs(...)'),
        # it never clears a pre-existing one. Pose ids are renumbered 1..N on every run,
        # so a partial relaunch could otherwise leave a stale model_merge_<id>.pdb paired
        # with a fresh score line for that same id in createOutputStep.
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
        """ Refines the filtered poses with RosettaDock, via run_protac_model.py
        --phase refine -> PROTAC-Model's own ros.rosetta(). Only inserted when doRefine
        is set (see _insertAllSteps). """
        # extra/rosetta/ must be a sibling of extra/frodock/: ros.rosetta() uses
        # hardcoded relative paths like '../frodock/...' to reach filterPosesStep's output.
        rosettaDir = self._getExtraPath('rosetta')
        os.makedirs(rosettaDir, exist_ok=True)
        # Same reasoning as filterPosesStep's cleanPath: ros.rosetta() only creates
        # extra/rosetta_results if missing, never clears a pre-existing one.
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
        """ Collects the final (filtered, optionally refined) ternary complex models into
        a SetOfAtomStructsChem, one AtomStruct per model, with its VoroMQA interface score
        (more negative = better) as a dynamic attribute. """
        if self.doRefine.get():
            resultsDir = self._getExtraPath('rosetta_results', 'all')
            resultsFile = os.path.join(resultsDir, 'results_rosetta.txt')
        else:
            resultsDir = self._getExtraPath('frodock_results', 'all')
            resultsFile = os.path.join(resultsDir, 'results_frodock.txt')

        # Without this check, a filtering/refinement phase that produced nothing (see the
        # two causes below) would silently fall through to an empty outputSet and a green
        # protocol - worse than a clear failure, since nothing downstream would flag it.
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
                # No header, 2 columns: "<pose_id> <score>" (verified against
                # PROTAC-Model's own utils/frodock.py and utils/rosetta.py).
                if not line.strip():
                    continue
                poseId, score = line.split()[:2]
                pdbFile = os.path.join(resultsDir, f'model_merge_{poseId}.pdb')
                if not os.path.exists(pdbFile):
                    self.info(f'Skipping pose {poseId}: {pdbFile} not found.')
                    continue
                atomStruct = AtomStruct(filename=pdbFile)
                atomStruct.setObjLabel(f'model_merge_{poseId}')
                # Dynamic attribute (AtomStruct has no built-in score field), same
                # mechanism protocol_flexDDG.py uses for its own per-item scores; name
                # matches pwchem's own '_score' convention (pwchem.objects.base).
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
        # Only when actually used: e3Ligand1/e3Ligand2 influence which poses survive
        # filtering (see _getLigLocateNum), so they belong in the provenance graph too.
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

        # PROTAC-Model's own generate_rosetta_para() (utils/rosetta.py) only assigns
        # gen_conf_num when cpu > 10, and then uses it unconditionally for -nstruct -
        # cpu<=10 is an UnboundLocalError, but only after frodockStep/filterPosesStep
        # have already run to completion. Caught here instead, before any of that runs.
        if self.doRefine.get() and self.numberOfThreads.get() <= 10:
            errors.append('"Refine with RosettaDock" needs more than 10 threads: '
                          'PROTAC-Model\'s own refinement code (generate_rosetta_para) '
                          'crashes with cpu<=10 (UnboundLocalError on gen_conf_num), and '
                          'that number also sets mpirun\'s -np and Rosetta\'s nstruct, so '
                          f'it cannot be silently raised for you. Got {self.numberOfThreads.get()} '
                          'threads.')

        # Both raise FileNotFoundError on the first missing tool home rather than a list,
        # so only one error is reported per call - acceptable, fixed one at a time.
        for check in (Plugin.getProtacModelEnviron, Plugin.getProtacModelPython):
            try:
                check()
            except FileNotFoundError as e:
                errors.append(str(e))

        return errors

    def _summary(self):
        summary = []
        # Set/Object.get() is for scalar attributes and always returns None on a Set -
        # isFinished() plus the hasattr check (set by _defineOutputs) is the right test,
        # matching protocol_flexDDG.py's own _summary().
        if self.isFinished() and hasattr(self, 'outputTernaryModels'):
            models = self.outputTernaryModels
            # createOutputStep now fails the protocol on 0 models (see its own comment),
            # so an empty outputTernaryModels shouldn't reach here - kept anyway, since a
            # finished protocol with no models used to crash _summary outright otherwise.
            if len(models):
                bestScore = min(model._score.get() for model in models)
                summary.append(f'Generated {len(models)} ternary complex model(s); best '
                               f'(most negative) interface score: {bestScore:.2f}.')
        return summary

    def _citations(self):
        return ['LeaverFay2011']

    # --------------------------- UTILS functions ------------------------------
    @staticmethod
    def _cleanReceptorOrTarget(inFile, outFile, ligandNameParam):
        """ Writes outFile as a water-stripped copy of inFile, additionally keeping only
        the named warhead among heteroatoms when ligandNameParam is set. Shared by
        convertInputStep for both the receptor and the target, which need the same
        cleanup logic. .upper() on the residue name since PDB resnames are uppercase and
        a hand-typed 'b96' otherwise wouldn't match CleanStructureSelect's het2keep. """
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
        """ Cleaned receptor PDB written by convertInputStep. Recomputed from extraPath
        (not cached as an instance attribute) since each step can run in its own
        process. """
        return self._getExtraPath('receptor.pdb')

    def _getTargetFile(self):
        """ Cleaned target PDB written by convertInputStep. Same recompute-not-cache
        reasoning as _getReceptorFile. """
        return self._getExtraPath('target.pdb')

    def _getFrodockBinDir(self):
        """ Directory holding the FRODOCK binary shim built by convertInputStep (see
        Plugin.prepareFrodockBinDir). Same recompute-not-cache reasoning as
        _getReceptorFile/_getTargetFile.
        Project-relative, like the two above: it is only ever consumed either in-process
        (prepareFrodockBinDir, which runs from the project directory) or as the FRODOCK
        environment variable, which Plugin.getProtacModelEnviron() makes absolute. """
        return self._getExtraPath('frodock_bin')

    def _getLigLocateNum(self):
        """ 2 when both E3 ligand conformers are given (ambiguous anchoring orientation,
        e.g. thalidomide-based degraders), 1 otherwise - matches PROTAC-Model's own
        main.py logic (lig_locate_num). """
        if self.e3Ligand1.get() is not None and self.e3Ligand2.get() is not None:
            return 2
        return 1

    @staticmethod
    def _getSmiArg(smilesParam):
        """ 'none' is PROTAC-Model's own convention (main.py/utils.frodock) for "no
        SMILES given" - it isn't a value a real SMILES string could take, so a stripped,
        non-empty param value is passed through unchanged. """
        value = smilesParam.get()
        return value.strip() if value else 'none'


