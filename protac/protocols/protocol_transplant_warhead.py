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
Copies a bound warhead from a holo structure onto a homologous apo structure by
superposing their binding pockets. The geometry lives in protac.utils.transplant.
"""

import dataclasses
import json

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
    Transplants a bound warhead from a source structure onto a homologous apo target,
    superposing only the pocket residues, matched by sequence alignment. Both pockets
    must be in an equivalent conformation; a mismatch shows up as clashes.
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
                       help='Homologous structure that receives the warhead.')
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
                            'equivalent in the target.')

    @classmethod
    def validateInstallation(cls):
        """ Needs no external tool, unlike the plugin-wide check. """
        return []

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Form values go as step arguments so that editing them reruns the steps.
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
        """ Keeps one chain (plus the warhead in the source) and rewrites modified
        residues as ATOM. They go in het2keep because cleanPDB treats them as HETATM. """
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
        """ Superposes, writes the ligand into the target and warns about clashes
        without failing. """
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
        mergePDBs(self._getFinalFile('target'), self._getTransplantedLigandFile(),
                 self._getOutputFile(), hetatm2=True)

        with open(self._getReportFile(), 'w') as f:
            json.dump(dataclasses.asdict(report), f, indent=2)

        if report.nClashes > 0:
            self.warning(
                f'{report.nClashes} clash(es) (< {clashCutoff} A) between the '
                f'transplanted warhead and the target structure (min distance '
                f'{report.minLigProtDistance:.2f} A). The target pocket is probably in '
                'a different conformation from the source one.')

    def createOutputStep(self):
        """ Output AtomStruct, with the warhead centroid stored as _siteCoords. """
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
