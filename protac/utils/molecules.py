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


def findMolByName(smallMolSet, name):
    """ The molecule called name in smallMolSet, cloned because sets reuse one instance
    while iterating. Fails when the name is missing or shared by several molecules (e.g.
    docking poses of one ligand) instead of silently picking one. """
    matches = [mol.clone() for mol in smallMolSet if mol.getMolName() == name]
    if not matches:
        raise ValueError(f'no molecule named "{name}" in that set.')
    if len(matches) > 1:
        raise ValueError(f'{len(matches)} molecules are named "{name}" in that set (for '
                         'instance several docking poses). Keep only the one to use, e.g. '
                         'with the "Operate set" protocol, and pick that set instead.')
    return matches[0]

