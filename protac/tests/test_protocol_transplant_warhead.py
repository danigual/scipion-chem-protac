# **************************************************************************
# *
# * Name:     test of protocol_transplant_warhead.py
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
Integration test for ProtPROTACTransplantWarhead: runs the actual Scipion protocol (not
just the pure protac.utils.transplant functions, already covered by
test_transplant_utils.py) against the same synthetic fixtures, inside a real Scipion
project. No external binaries/*_HOME variables are needed (everything downstream of
pwchem.utils.cleanPDB/mergePDBs is Biopython/numpy), so this never skips.
"""

import os

from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import ProtImportPdb

from protac.protocols.protocol_transplant_warhead import ProtPROTACTransplantWarhead

_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'transplant')
_SOURCE_PDB = os.path.join(_DATA_DIR, 'toy_source.pdb')
_TARGET_PDB = os.path.join(_DATA_DIR, 'toy_target.pdb')


class TestTransplantWarhead(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        protImportSource = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                           pdbFile=_SOURCE_PDB)
        cls.launchProtocol(protImportSource)
        cls.sourceStructure = protImportSource.outputPdb

        protImportTarget = cls.newProtocol(ProtImportPdb, inputPdbData=1,
                                           pdbFile=_TARGET_PDB)
        cls.launchProtocol(protImportTarget)
        cls.targetStructure = protImportTarget.outputPdb

    def test_happy_path(self):
        args = {'sourceStructure': self.sourceStructure, 'sourceChain': 'A',
               'sourceLigandName': 'LIG', 'targetStructure': self.targetStructure,
               'targetChain': 'A', 'minPocketPairs': 4}
        prot = self.newProtocol(ProtPROTACTransplantWarhead, **args)
        self.launchProtocol(prot)

        outputStructure = getattr(prot, 'outputStructure', None)
        self.assertIsNotNone(outputStructure, 'No outputStructure was generated')
        self.assertTrue(os.path.exists(outputStructure.getFileName()))

        siteCoords = outputStructure._siteCoords.get()
        parts = siteCoords.split(',')
        self.assertEqual(len(parts), 3)
        for part in parts:
            float(part)  # raises if not a valid coordinate

        # The warhead is written into the target, on the target chain.
        with open(outputStructure.getFileName()) as f:
            ligLines = [l for l in f if l.startswith('HETATM') and l[17:20] == 'LIG']
        self.assertEqual(len(ligLines), 4)
        self.assertTrue(all(l[21] == 'A' for l in ligLines))

    def test_min_pocket_pairs_too_high_fails(self):
        args = {'sourceStructure': self.sourceStructure, 'sourceChain': 'A',
               'sourceLigandName': 'LIG', 'targetStructure': self.targetStructure,
               'targetChain': 'A', 'minPocketPairs': 100}
        prot = self.newProtocol(ProtPROTACTransplantWarhead, **args)
        with self.assertRaises(Exception):
            self.launchProtocol(prot)
        # Fails for the right reason, not just anywhere.
        self.assertTrue(prot.isFailed())
        self.assertIn('need >= 100', prot.getErrorMessage())
