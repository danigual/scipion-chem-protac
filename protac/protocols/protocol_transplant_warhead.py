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
Transplants a small-molecule warhead from a "source" structure that already has it
bound onto a homologous "target" structure that doesn't (e.g. an apo crystal structure
of the same protein, or of a different isoform), by superposing the two on their binding
pocket - restricted to the residues near the warhead in the source, matched to the
target by sequence alignment rather than residue numbering, since homologs commonly
have different residue numbering.

Thin wrapper around protac.utils.transplant, which does all the actual geometry/PDB
text work: this module only translates Scipion inputs (AtomStruct pointers, form
parameters) into file paths and the module's TransplantReport into Scipion objects.
"""

import json
import os

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.utils import cleanPDB, mergePDBs

from protac.constants import (POCKET_CUTOFF, CLASH_CUTOFF, MIN_POCKET_PAIRS,
                              MODRES_TO_CANONICAL)
from protac.utils.transplant import (runTransplant, rewriteModifiedResiduesAsAtom,
                                     writeTransplantedLigand)


class ProtPROTACTransplantWarhead(EMProtocol):
    """
    Transplants a bound warhead from a source structure onto a homologous, apo target
    structure, via sequence-alignment-guided, pocket-restricted structural superposition.

    Scope: only transplant between structurally equivalent conformations. A genuine
    conformational mismatch between the warhead and the target pocket (e.g. a type-II
    inhibitor that needs a DFG-out-like conformation, transplanted onto a DFG-in target)
    is expected to surface as clashes in the quality control below, but this protocol
    does not attempt to diagnose or resolve the mismatch itself.
    """
    _label = 'Warhead transplant (homology-guided)'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        group = form.addGroup('Source (holo)')
        group.addParam('sourceStructure', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Source structure (holo)',
                       help='Structure already bound to the warhead to transplant.')
        group.addParam('sourceChain', params.StringParam, default='A', allowsNull=False,
                       label='Source chain',
                       help='Chain ID of the protein in the source structure, as it '
                            'appears in the PDB.')
        group.addParam('sourceLigandName', params.StringParam, allowsNull=False,
                       label='Warhead residue name (source)',
                       help='Residue name of the bound warhead in the source structure, '
                            'as it appears in the PDB HETATM records (e.g. "B96").')

        group = form.addGroup('Target (apo)')
        group.addParam('targetStructure', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=False, label='Target structure (apo)',
                       help='Homologous structure to receive the transplanted warhead. '
                            'Its pocket residues (matched to the source pocket by '
                            'sequence alignment, not by residue number) are used to '
                            'restrict the superposition.')
        group.addParam('targetChain', params.StringParam, default='A', allowsNull=False,
                       label='Target chain',
                       help='Chain ID of the protein in the target structure, as it '
                            'appears in the PDB.')

        group = form.addGroup('Quality control')
        group.addParam('pocketCutoff', params.FloatParam, default=POCKET_CUTOFF,
                       expertLevel=params.LEVEL_ADVANCED, label='Pocket cutoff (A)',
                       help='A source residue counts as pocket if it has an atom within '
                            'this distance of the warhead.')
        group.addParam('clashCutoff', params.FloatParam, default=CLASH_CUTOFF,
                       expertLevel=params.LEVEL_ADVANCED, label='Clash cutoff (A)',
                       help='Ligand-protein contacts closer than this, in the target, '
                            'are counted as clashes.')
        group.addParam('minPocketPairs', params.IntParam, default=MIN_POCKET_PAIRS,
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Minimum matched pocket residues',
                       help='Aborts if fewer source pocket residues have an aligned '
                            'equivalent in the target: the superposition would not be '
                            'reliable below this.')

    @classmethod
    def validateInstallation(cls):
        """ Overrides the plugin-wide default, which checks Plugin._pathVars -
        PROTAC_MODEL_HOME, a dependency of ProtPROTACModel that this protocol never
        touches. Everything here (pwchem.utils, Biopython,
        numpy) is already part of the scipion3 environment - nothing to validate. """
        return []

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Every value that changes prepareInputsStep/transplantStep's output is passed
        # as a funcArgs so Scipion's funcName+argsStr equality check (protocol.py) can
        # tell a re-run with edited form values apart from a genuine resume - otherwise
        # "Continue" after fixing e.g. a typo'd sourceLigandName would silently reuse the
        # stale extra/ files from the first, wrong run.
        prepArgs = [self.sourceChain.get(), self.targetChain.get(),
                   self.sourceLigandName.get()]
        prepId = self._insertFunctionStep(self.prepareInputsStep, *prepArgs,
                                          prerequisites=[])

        transplantArgs = prepArgs + [self.pocketCutoff.get(), self.clashCutoff.get(),
                                     self.minPocketPairs.get()]
        transplantId = self._insertFunctionStep(self.transplantStep, *transplantArgs,
                                                prerequisites=[prepId])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[transplantId])

    def prepareInputsStep(self, sourceChain, targetChain, sourceLigandName):
        """ Cleans the source/target PDBs down to a single protein chain (+ the named
        warhead, for the source), then rewrites modified residues (PTR/TPO/SEP/MSE) as
        ATOM. Modified residues are kept in cleanPDB's het2keep on purpose: cleanPDB
        treats them as heteroatoms too (they aren't among the 20 standard residues
        Biopython recognises), so without this they would be stripped out before
        rewriteModifiedResiduesAsAtom ever saw them. """
        sourceChain, targetChain = sourceChain.strip(), targetChain.strip()
        modResNames = list(MODRES_TO_CANONICAL)
        cleanPDB(self.sourceStructure.get().getFileName(), self._getCleanedFile('source'),
                 waters=True, hetatm=True,
                 het2keep=[sourceLigandName.strip().upper()] + modResNames,
                 chainIds=[sourceChain])
        cleanPDB(self.targetStructure.get().getFileName(), self._getCleanedFile('target'),
                 waters=True, hetatm=True, het2keep=modResNames,
                 chainIds=[targetChain])

        rewriteModifiedResiduesAsAtom(self._getCleanedFile('source'),
                                      self._getFinalFile('source'), sourceChain)
        rewriteModifiedResiduesAsAtom(self._getCleanedFile('target'),
                                      self._getFinalFile('target'), targetChain)

    def transplantStep(self, sourceChain, targetChain, sourceLigandName, pocketCutoff,
                       clashCutoff, minPocketPairs):
        """ Runs the pure alignment/superposition/QC algorithm, writes the transplanted
        ligand and merges it into the target, and reports any clash - without aborting
        the protocol, since a clash is a QC signal for the user to act on, not by
        itself a broken run (matches the reference approach this generalizes). """
        sourceChain, targetChain = sourceChain.strip(), targetChain.strip()
        sourceLigandName = sourceLigandName.strip().upper()
        newCoords, report = runTransplant(
            self._getFinalFile('source'), sourceChain, sourceLigandName,
            self._getFinalFile('target'), targetChain,
            pocketCutoff=pocketCutoff, clashCutoff=clashCutoff,
            minPocketPairs=minPocketPairs)

        writeTransplantedLigand(self._getFinalFile('source'), sourceChain,
                                sourceLigandName, targetChain,
                                newCoords, self._getTransplantedLigandFile())
        # hetatm2=True is defensive/documentational here, not load-bearing:
        # writeTransplantedLigand only ever writes HETATM lines, and mergePDBs only
        # rewrites lines starting with 'ATOM' - it is a no-op on this input, kept so the
        # call reads correctly if writeTransplantedLigand's output format ever changes.
        mergePDBs(self._getFinalFile('target'), self._getTransplantedLigandFile(),
                 self._getOutputFile(), hetatm2=True)

        with open(self._getReportFile(), 'w') as f:
            json.dump(report.asdict(), f, indent=2)

        if report.nClashes > 0:
            self.warning(
                f'{report.nClashes} clash(es) (< {clashCutoff} A) between the '
                f'transplanted warhead and the target structure (min distance '
                f'{report.minLigProtDistance:.2f} A). This structure likely fails the '
                'acceptance criterion for downstream use. A frequent cause is a '
                'conformational mismatch between the warhead and the target pocket '
                "(e.g. a type-II inhibitor that needs a DFG-out-like conformation, "
                'transplanted onto a DFG-in target) - this check flags the symptom '
                'generically; it does not diagnose the specific structural motif '
                'involved.')

    def createOutputStep(self):
        """ Packages the merged PDB as an AtomStruct, with the transplanted warhead's
        centroid as a dynamic '_siteCoords' attribute so a future wizard
        can fill in ProtPROTACModel's manual siteCoords field
        from it without re-parsing the structure. """
        with open(self._getReportFile()) as f:
            report = json.load(f)

        atomStruct = AtomStruct(filename=self._getOutputFile())
        x, y, z = report['centroid']
        atomStruct._siteCoords = pwobj.String(f'{x:.3f},{y:.3f},{z:.3f}')

        self._defineOutputs(outputStructure=atomStruct)
        self._defineSourceRelation(self.sourceStructure, atomStruct)
        self._defineSourceRelation(self.targetStructure, atomStruct)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        if not self.sourceChain.get().strip():
            errors.append('"Source chain" cannot be empty.')
        if not self.targetChain.get().strip():
            errors.append('"Target chain" cannot be empty.')
        if not self.sourceLigandName.get().strip():
            errors.append('"Warhead residue name (source)" cannot be empty: there is '
                          'nothing to transplant without it.')
        if self.pocketCutoff.get() <= 0:
            errors.append('"Pocket cutoff (A)" must be positive.')
        if self.clashCutoff.get() <= 0:
            errors.append('"Clash cutoff (A)" must be positive.')
        if self.minPocketPairs.get() < 1:
            errors.append('"Minimum matched pocket residues" must be at least 1.')
        return errors

    def _summary(self):
        summary = []
        if self.isFinished() and hasattr(self, 'outputStructure'):
            with open(self._getReportFile()) as f:
                report = json.load(f)
            summary.append(
                f'Sequence identity: {report["identityPct"]:.1f}%. Pocket: '
                f'{report["nPocketPairs"]}/{report["nPocketSource"]} residue(s) matched, '
                f'RMSD {report["pocketRmsd"]:.2f} A. Clashes: {report["nClashes"]} '
                f'(min distance {report["minLigProtDistance"]:.2f} A).')
        return summary

    def _citations(self):
        # TODO: protac/bibtex.py doesn't exist yet, so any key
        # returned here is silently dropped regardless - leaving this empty rather than
        # guessing an unverified Biopython citation key.
        return []

    # --------------------------- UTILS functions ------------------------------
    def _getCleanedFile(self, which):
        return self._getExtraPath(f'{which}_clean.pdb')

    def _getFinalFile(self, which):
        return self._getExtraPath(f'{which}_final.pdb')

    def _getTransplantedLigandFile(self):
        return self._getExtraPath('ligand_transplanted.pdb')

    def _getOutputFile(self):
        return self._getExtraPath('transplanted.pdb')

    def _getReportFile(self):
        return self._getExtraPath('report.json')
