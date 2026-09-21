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
Unit tests for protac.utils.transplant: pure functions (Biopython/numpy only), no
Scipion project needed, so these never skip - unlike protac's other external-tool
dependencies, nothing here requires a binary or a *_HOME environment variable.

Fixtures (protac/tests/data/transplant/toy_source.pdb, toy_target.pdb) are synthetic,
never real p38/CRBN data: a straight 15-residue toy chain with a 4-atom "LIG" ligand
near residues 6-9 and one HETATM PTR at residue 3 (toy_source.pdb), and the same chain
rigidly transformed by a known rotation/translation with one conservative mutation at
residue 1, no ligand (toy_target.pdb). The known transform lets test_runTransplant_
happy_path check the recovered centroid against an exact expected value instead of just
"looks reasonable".
"""

import os
import tempfile
import unittest

import numpy as np

from protac.utils import transplant as T

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'transplant')
_SOURCE_PDB = os.path.join(_DATA_DIR, 'toy_source.pdb')
_TARGET_PDB = os.path.join(_DATA_DIR, 'toy_target.pdb')

# Exact rotation/translation used to build toy_target.pdb from toy_source.pdb's backbone
# (see the generation script used when the fixtures were created) - reused here to check
# runTransplant's output against a known value, not just plausibility.
_THETA = np.radians(30)
_ROT = np.array([[np.cos(_THETA), -np.sin(_THETA), 0],
                 [np.sin(_THETA), np.cos(_THETA), 0],
                 [0, 0, 1]])
_TRAN = np.array([10.0, 5.0, 2.0])


def _parsedStructures():
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    return (parser.get_structure('src', _SOURCE_PDB),
           parser.get_structure('tgt', _TARGET_PDB))


class TestTransplantUtils(unittest.TestCase):
    def test_chainResidues_excludes_ligand_keeps_modified_residue(self):
        src, _ = _parsedStructures()
        residues = T.chainResidues(src, 'A')
        resnames = [r.get_resname() for r in residues]
        self.assertEqual(len(residues), 15)
        self.assertIn('PTR', resnames)
        self.assertNotIn('LIG', resnames)

    def test_alignedPairs_identity_and_pairs(self):
        src, tgt = _parsedStructures()
        srcRes = T.chainResidues(src, 'A')
        tgtRes = T.chainResidues(tgt, 'A')
        aln = T.alignedPairs(srcRes, tgtRes)
        # No gaps expected (same length, no indels): all 15 positions pair up.
        self.assertEqual(len(aln.pairs), 15)
        # 2 of 15 positions differ: the explicit mutation at residue 1, and residue 3
        # (PTR in the source vs its unmodified residue in the target) - 13/15.
        self.assertAlmostEqual(aln.identityPct, 100.0 * 13 / 15, places=1)

    def test_pocketResidues_cutoff(self):
        src, _ = _parsedStructures()
        srcRes = T.chainResidues(src, 'A')
        ligAtoms = T.ligandAtoms(src, 'A', 'LIG')
        pocket = T.pocketResidues(srcRes, ligAtoms, cutoff=6.0)
        self.assertEqual(len(pocket), 6)
        self.assertEqual(T.pocketResidues(srcRes, ligAtoms, cutoff=0.01), set())

    def test_superposePocket_raises_below_min_pairs(self):
        src, tgt = _parsedStructures()
        srcRes = T.chainResidues(src, 'A')
        tgtRes = T.chainResidues(tgt, 'A')
        ligAtoms = T.ligandAtoms(src, 'A', 'LIG')
        aln = T.alignedPairs(srcRes, tgtRes)
        pocket = T.pocketResidues(srcRes, ligAtoms)
        with self.assertRaises(ValueError):
            T.superposePocket(aln.pairs, pocket, minPairs=100)

    def test_superposePocket_rmsd_near_zero(self):
        src, tgt = _parsedStructures()
        srcRes = T.chainResidues(src, 'A')
        tgtRes = T.chainResidues(tgt, 'A')
        ligAtoms = T.ligandAtoms(src, 'A', 'LIG')
        aln = T.alignedPairs(srcRes, tgtRes)
        pocket = T.pocketResidues(srcRes, ligAtoms)
        sup = T.superposePocket(aln.pairs, pocket, minPairs=4)
        # None of the pocket residues (near residues 6-11) is the mutated residue 1, so
        # the pocket is a pure rigid-body copy: RMSD should be ~0 (PDB text rounds
        # coordinates to 3 decimals, hence the small but nonzero tolerance).
        self.assertLess(sup.rms, 0.01)

    def test_checkClashes_detects_and_counts(self):
        newLigCoords = np.array([[0.0, 0.0, 0.0]])
        overlapping = np.array([[0.1, 0.0, 0.0]])
        farAway = np.array([[100.0, 0.0, 0.0]])
        self.assertEqual(T.checkClashes(newLigCoords, overlapping, clashCutoff=2.0).nClashes, 1)
        self.assertEqual(T.checkClashes(newLigCoords, farAway, clashCutoff=2.0).nClashes, 0)

    def test_runTransplant_happy_path(self):
        newCoords, report = T.runTransplant(_SOURCE_PDB, 'A', 'LIG', _TARGET_PDB, 'A',
                                            minPocketPairs=4)
        self.assertEqual(report.nClashes, 0)
        self.assertLess(report.pocketRmsd, 0.01)

        srcStruct, _ = _parsedStructures()
        origLigCoords = np.array([a.coord for a in T.ligandAtoms(srcStruct, 'A', 'LIG')])
        expectedCentroid = origLigCoords.mean(axis=0).dot(_ROT.T) + _TRAN
        np.testing.assert_allclose(report.centroid, expectedCentroid, atol=0.01)
        np.testing.assert_allclose(newCoords.mean(axis=0), expectedCentroid, atol=0.01)

    def test_runTransplant_no_ligand_atoms_raises(self):
        with self.assertRaises(ValueError):
            T.runTransplant(_SOURCE_PDB, 'A', 'NOPE', _TARGET_PDB, 'A')

    def test_rewriteModifiedResiduesAsAtom_rewrites_and_filters(self):
        with tempfile.TemporaryDirectory() as tmp:
            outPdb = os.path.join(tmp, 'rewritten.pdb')
            T.rewriteModifiedResiduesAsAtom(_SOURCE_PDB, outPdb, 'A')
            with open(outPdb) as f:
                lines = [l for l in f if l.startswith(('ATOM', 'HETATM'))]
            ptrLines = [l for l in lines if l[17:20].strip() == 'PTR']
            self.assertEqual(len(ptrLines), 1)
            self.assertTrue(ptrLines[0].startswith('ATOM'))
            # The ligand is untouched (still HETATM, not a modified residue).
            ligLines = [l for l in lines if l[17:20].strip() == 'LIG']
            self.assertTrue(all(l.startswith('HETATM') for l in ligLines))

    def test_rewriteModifiedResiduesAsAtom_filters_altloc_and_hydrogen(self):
        # toy_source.pdb has no altlocs/hydrogens, so this targeted edge case is built
        # inline instead of adding another fixture file just for it.
        content = (
            'ATOM      1  CA AALA A   1      10.000  10.000  10.000  0.50  0.00           C\n'
            'ATOM      2  CA BALA A   1      11.000  10.000  10.000  0.50  0.00           C\n'
            'ATOM      3  H   ALA A   1      10.500  10.500  10.000  1.00  0.00           H\n'
            'END\n'
        )
        with tempfile.TemporaryDirectory() as tmp:
            inPdb = os.path.join(tmp, 'in.pdb')
            outPdb = os.path.join(tmp, 'out.pdb')
            with open(inPdb, 'w') as f:
                f.write(content)
            T.rewriteModifiedResiduesAsAtom(inPdb, outPdb, 'A')
            with open(outPdb) as f:
                lines = [l for l in f if l.startswith(('ATOM', 'HETATM'))]
            # Only the primary altloc ('A') survives, and the hydrogen is dropped.
            self.assertEqual(len(lines), 1)
            self.assertEqual(lines[0][76:78].strip(), 'C')

    def test_writeTransplantedLigand_reassigns_chain_and_coords(self):
        srcStruct, _ = _parsedStructures()
        ligAtoms = T.ligandAtoms(srcStruct, 'A', 'LIG')
        newCoords = np.array([[float(i), 0.0, 0.0] for i in range(len(ligAtoms))])
        with tempfile.TemporaryDirectory() as tmp:
            outPdb = os.path.join(tmp, 'lig.pdb')
            T.writeTransplantedLigand(_SOURCE_PDB, 'A', 'LIG', 'B', newCoords, outPdb)
            with open(outPdb) as f:
                lines = [l for l in f if l.startswith('HETATM')]
            self.assertEqual(len(lines), len(ligAtoms))
            for i, line in enumerate(lines):
                self.assertEqual(line[21], 'B')
                self.assertAlmostEqual(float(line[30:38]), float(i), places=2)

            with self.assertRaises(ValueError):
                T.writeTransplantedLigand(_SOURCE_PDB, 'A', 'LIG', 'B', newCoords[:-1], outPdb)


if __name__ == '__main__':
    unittest.main()
