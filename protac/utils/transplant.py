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
Pure Python (no pyworkflow/pwem, only pwchem.utils for PDB line helpers and the
standard-residue table) implementation of homology-guided warhead transplant: given a
"source" structure with a bound small-molecule warhead and a homologous, apo "target"
structure, superpose the two on their binding pocket (residues near the warhead in the
source, matched to the target by sequence alignment rather than residue numbering) and
carry the warhead's coordinates over.

Kept deliberately outside protac/protocols/: everything here is plain geometry/text
manipulation, testable without a Scipion project. The protocol that wraps this
(ProtPROTACTransplantWarhead) is a thin layer translating Scipion inputs to file paths
and this module's outputs to Scipion objects.

Scope: only "transplant between equivalent conformations". A genuine conformational
mismatch (e.g. a type-II warhead that needs a DFG-out-like pocket, transplanted onto a
DFG-in target) is expected to surface as clashes here, but diagnosing *why* is out of
scope - see ProtPROTACTransplantWarhead's clash warning.
"""

from dataclasses import dataclass

import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.Align import PairwiseAligner, substitution_matrices

from pwchem.utils import RESIDUES3TO1, splitPDBLine, writePDBLine

from protac.constants import (MODRES_TO_CANONICAL, POCKET_CUTOFF, CLASH_CUTOFF,
                              MIN_POCKET_PAIRS)

# Residues that count as "main chain" for alignment purposes: the 20 standard amino
# acids plus the handful of modified residues (PTR/TPO/SEP/MSE) that PDB files often
# record as HETATM even though they are part of the sequence.
_CHAIN_RESNAMES = {**RESIDUES3TO1, **MODRES_TO_CANONICAL}


def _resLetter(resname):
    """ One-letter alignment code (canonical, uppercase) for a chain residue name. """
    return _CHAIN_RESNAMES[resname]


def chainResidues(structure, chainId):
    """ Amino-acid residues (standard or modified, see _CHAIN_RESNAMES) of a chain, in
    order, restricted to those with a CA atom. """
    chain = structure[0][chainId]
    return [res for res in chain if res.get_resname() in _CHAIN_RESNAMES and 'CA' in res]


def ligandAtoms(structure, chainId, resname):
    """ Heavy (non-hydrogen) atoms of the ligand named `resname` in `chainId`. """
    atoms = []
    for res in structure[0][chainId]:
        if res.get_resname() == resname:
            atoms.extend(a for a in res if a.element != 'H')
    return atoms


@dataclass
class AlignmentResult:
    pairs: list           # [(source_residue, target_residue), ...], gaps excluded
    identityPct: float


def alignedPairs(srcRes, tgtRes):
    """ Global BLOSUM62 alignment of the two residue lists' sequences, returned as
    pairs of aligned Bio.PDB.Residue objects (positions aligned to a gap are dropped). """
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = 'global'

    srcSeq = ''.join(_resLetter(r.get_resname()) for r in srcRes)
    tgtSeq = ''.join(_resLetter(r.get_resname()) for r in tgtRes)
    aln = aligner.align(srcSeq, tgtSeq)[0]

    pairs = []
    for (s0, s1), (t0, t1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(s1 - s0):
            pairs.append((srcRes[s0 + k], tgtRes[t0 + k]))

    identity = sum(1 for a, b in pairs if a.get_resname() == b.get_resname())
    identityPct = 100.0 * identity / max(len(srcSeq), len(tgtSeq)) if (srcSeq or tgtSeq) else 0.0
    return AlignmentResult(pairs=pairs, identityPct=identityPct)


def pocketResidues(srcRes, ligAtoms, cutoff=POCKET_CUTOFF):
    """ Subset of srcRes with at least one heavy atom within `cutoff` A of ligAtoms. """
    ligCoords = np.array([a.coord for a in ligAtoms])
    pocket = set()
    for res in srcRes:
        coords = np.array([a.coord for a in res if a.element != 'H'])
        d = np.linalg.norm(coords[:, None, :] - ligCoords[None, :, :], axis=-1)
        if d.min() <= cutoff:
            pocket.add(res)
    return pocket


@dataclass
class SuperpositionResult:
    rotran: tuple          # (rot: np.ndarray[3,3], tran: np.ndarray[3]), source->target
    rms: float              # CA RMSD of the pocket after superposition
    nPocketPairs: int


def superposePocket(pairs, pocket, minPairs=MIN_POCKET_PAIRS):
    """ Restricts `pairs` to the ones whose source residue is in `pocket`, and
    superposes their CA atoms (target <- source). Raises ValueError if fewer than
    `minPairs` pocket residues have an aligned target equivalent - the superposition
    would not be reliable below that. """
    pocketPairs = [(a, b) for a, b in pairs if a in pocket]
    if len(pocketPairs) < minPairs:
        raise ValueError(
            f'Only {len(pocketPairs)} pocket residue(s) have an aligned equivalent in '
            f'the target (need >= {minPairs}): the superposition would not be reliable.')

    sup = Superimposer()
    sup.set_atoms([b['CA'] for _, b in pocketPairs], [a['CA'] for a, _ in pocketPairs])
    return SuperpositionResult(rotran=sup.rotran, rms=sup.rms, nPocketPairs=len(pocketPairs))


def transplantLigand(ligAtoms, rotran):
    """ Applies (rot, tran) to the coordinates of ligAtoms. Returns an [n, 3] array. """
    rot, tran = rotran
    coords = np.array([a.coord for a in ligAtoms])
    return np.dot(coords, rot) + tran


@dataclass
class ClashReport:
    minDistance: float
    nClashes: int


def checkClashes(newLigCoords, targetProtCoords, clashCutoff=CLASH_CUTOFF):
    """ Minimum ligand-protein distance and number of contacts below `clashCutoff` A,
    between the transplanted ligand and the target structure. """
    d = np.linalg.norm(np.array(newLigCoords)[:, None, :]
                       - np.array(targetProtCoords)[None, :, :], axis=-1)
    return ClashReport(minDistance=float(d.min()), nClashes=int((d < clashCutoff).sum()))


def centroid(coords):
    """ (x, y, z) centroid of an array/list of coordinates. """
    return tuple(np.asarray(coords).mean(axis=0))


@dataclass
class TransplantReport:
    """ Everything the algorithm can report about the run, structured so the protocol
    can log/summarize it without re-running or re-parsing anything. """
    identityPct: float
    nAlignedPairs: int
    nPocketSource: int
    nPocketPairs: int
    pocketRmsd: float
    minLigProtDistance: float
    nClashes: int
    centroid: tuple



def runTransplant(sourcePdb, sourceChain, ligResname, targetPdb, targetChain,
                  pocketCutoff=POCKET_CUTOFF, clashCutoff=CLASH_CUTOFF,
                  minPocketPairs=MIN_POCKET_PAIRS):
    """ Runs the full pure-geometry pipeline (parse -> align -> pocket -> superpose ->
    transplant -> QC). Returns (newLigCoords: np.ndarray[n,3], report: TransplantReport).

    Does not write any output file: building the final merged PDB mixes these new
    coordinates with text already written elsewhere (chain reassignment, HETATM
    formatting) - see writeTransplantedLigand() and pwchem.utils.mergePDBs, both meant
    to be called by the protocol after this. """
    parser = PDBParser(QUIET=True)
    src = parser.get_structure('src', sourcePdb)
    tgt = parser.get_structure('tgt', targetPdb)

    srcRes = chainResidues(src, sourceChain)
    tgtRes = chainResidues(tgt, targetChain)
    ligAtoms = ligandAtoms(src, sourceChain, ligResname)
    if not ligAtoms:
        raise ValueError(
            f'No {ligResname} HETATM records found in chain {sourceChain} of the source '
            'structure.')

    aln = alignedPairs(srcRes, tgtRes)
    pocket = pocketResidues(srcRes, ligAtoms, cutoff=pocketCutoff)
    sup = superposePocket(aln.pairs, pocket, minPairs=minPocketPairs)

    newLigCoords = transplantLigand(ligAtoms, sup.rotran)
    tgtProtCoords = np.array([a.coord for r in tgtRes for a in r if a.element != 'H'])
    clashes = checkClashes(newLigCoords, tgtProtCoords, clashCutoff=clashCutoff)

    report = TransplantReport(
        identityPct=aln.identityPct,
        nAlignedPairs=len(aln.pairs),
        nPocketSource=len(pocket),
        nPocketPairs=sup.nPocketPairs,
        pocketRmsd=sup.rms,
        minLigProtDistance=clashes.minDistance,
        nClashes=clashes.nClashes,
        centroid=centroid(newLigCoords),
    )
    return newLigCoords, report


def rewriteModifiedResiduesAsAtom(inPdb, outPdb, chainId):
    """ Copies the ATOM/HETATM records of `chainId` from inPdb to outPdb, rewriting
    HETATM records of modified residues (MODRES_TO_CANONICAL) as ATOM - they are main
    chain, not ligands, and pwchem.utils.cleanPDB has no option to do this itself, so
    left as HETATM they would be dropped by any later hetatm=True cleaning and would
    otherwise contaminate a HETATM-only warhead centroid/selection. Also drops alternate
    (non-primary) altloc records and hydrogens: cleanPDB does not handle either.
    Any non-ATOM/HETATM line (headers, TER, END...) is passed through unchanged. """
    with open(inPdb) as fIn, open(outPdb, 'w') as fOut:
        for line in fIn:
            if not line.startswith(('ATOM', 'HETATM')):
                fOut.write(line)
                continue

            parts = splitPDBLine(line)
            recordType, _, _, resName, chain = parts[:5]
            elementSymbol = parts[-1]
            altLoc = line[16]

            if chain != chainId or altLoc not in (' ', 'A') or elementSymbol == 'H':
                continue

            if recordType == 'HETATM' and resName in MODRES_TO_CANONICAL:
                parts[0] = 'ATOM'
                line = writePDBLine(parts)
            fOut.write(line)


def writeTransplantedLigand(sourcePdb, sourceChain, ligResname, targetChain, newCoords,
                            outPdb):
    """ Writes the ligand's HETATM records (as read from sourcePdb, chain sourceChain,
    residue name ligResname - same selection criteria as ligandAtoms(), so the atom
    order matches `newCoords`) with their coordinates replaced by `newCoords` and their
    chain reassigned to `targetChain`, ready to be appended to the target structure via
    pwchem.utils.mergePDBs(..., hetatm2=True). """
    ligParts = []
    with open(sourcePdb) as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            parts = splitPDBLine(line)
            if parts[3] != ligResname or parts[4] != sourceChain:
                continue
            if line[16] not in (' ', 'A') or parts[-1] == 'H':
                continue
            ligParts.append(parts)

    if len(ligParts) != len(newCoords):
        raise ValueError(
            f'{len(ligParts)} {ligResname} atom(s) read from {sourcePdb}, but '
            f'{len(newCoords)} transplanted coordinate(s) were given - the source PDB '
            'changed between ligandAtoms() and writeTransplantedLigand().')

    with open(outPdb, 'w') as fOut:
        for parts, xyz in zip(ligParts, newCoords):
            parts = list(parts)
            parts[4] = targetChain
            parts[6], parts[7], parts[8] = f'{xyz[0]:.3f}', f'{xyz[1]:.3f}', f'{xyz[2]:.3f}'
            fOut.write(writePDBLine(parts))
