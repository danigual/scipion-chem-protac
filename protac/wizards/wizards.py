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

from pyworkflow.gui import dialog
from pwem.wizards import VariableWizard

from pwchem.wizards import SelectElementWizard, SelectLigandWizard

from protac.protocols import ProtPROTACModel, ProtPRosettaC
from protac.utils.molecules import ligandCentroid


class LigandCentroidWizard(VariableWizard):
    """ Fills a coordinates field with the centroid of a named ligand in a structure. """
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, outputParam = self.getInputOutput(form)
        struct = getattr(protocol, inputParam[0]).get()
        ligName = getattr(protocol, inputParam[1]).get()
        if struct is None or not ligName:
            dialog.showError('Missing input', 'Choose the structure and its warhead '
                             'residue name first.', form.root)
            return
        try:
            coords = ligandCentroid(struct.getFileName(), ligName)
        except ValueError as e:
            dialog.showError('Cannot compute the centroid', str(e), form.root)
            return
        form.setVar(outputParam[0], '%.2f,%.2f,%.2f' % coords)


# Warhead residue names, picked from the heteroatoms of each structure.
for structParam, nameParam in (('inputReceptor', 'receptorLigandName'),
                               ('inputTarget', 'targetLigandName')):
    SelectLigandWizard().addTarget(protocol=ProtPROTACModel, targets=[nameParam],
                                   inputs=[structParam], outputs=[nameParam])

LigandCentroidWizard().addTarget(protocol=ProtPROTACModel, targets=['siteCoords'],
                                 inputs=['inputReceptor', 'receptorLigandName'],
                                 outputs=['siteCoords'])

# Molecule names within the input sets.
for protocol, nameParam, namesMethod in (
        (ProtPROTACModel, 'e3Ligand1Name', 'getE3LigandNames'),
        (ProtPROTACModel, 'e3Ligand2Name', 'getE3LigandNames'),
        (ProtPRosettaC, 'head1Name', 'getHead1Names'),
        (ProtPRosettaC, 'head2Name', 'getHead2Names')):
    SelectElementWizard().addTarget(protocol=protocol, targets=[nameParam],
                                    inputs=[namesMethod], outputs=[nameParam])
