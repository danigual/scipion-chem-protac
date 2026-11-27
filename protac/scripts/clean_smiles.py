# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Gutiérrez (Internship student, github user: danigual)
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
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

'''Script to normalize the SMILES derived from ligand coordinate files, run in the
rdkit-env.

SMILES written from a molecule parsed out of a coordinate file often carry fixed
valences: OpenBabel fills in the valence field of the SDF atom block, RDKit reads it as
"this atom takes no implicit hydrogens" and writes every atom bracketed ([C], [N], [O]).
The skeleton is right but the molecule has no hydrogens at all, which is useless as a
bond order template: PROTAC-Model passes these SMILES to
AllChem.AssignBondOrdersFromTemplate, which then fails with "No matching found".

This script rewrites each SMILES of a csv (as written by rdkit_IO.py with
-of smiles_csv) letting RDKit work the hydrogen counts out from the perceived bonds,
and optionally checks that the warhead SMILES are substructures of the PROTAC, which
is what PROTAC-Model needs to grow the linker conformations.
'''

import sys, os, argparse, csv

from rdkit import Chem
from rdkit import RDLogger


def normalizeSmiles(smi):
    '''Rewrite a SMILES, dropping the fixed valences that come from coordinate files'''
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        raise ValueError('RDKit cannot parse the SMILES: {}'.format(smi))

    for atom in mol.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumExplicitHs(0)

    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol)

    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def readSmilesCsv(inFile):
    '''Parse the csv written by rdkit_IO.py, returning its rows and the smiles column'''
    with open(inFile, newline='') as f:
        reader = csv.DictReader(f)
        fieldNames = reader.fieldnames or []
        smilesColumn = next((col for col in fieldNames if 'smiles' in col.lower()), None)

        if smilesColumn is None:
            raise ValueError('No SMILES column found in {}. Columns: {}'.
                             format(inFile, fieldNames))

        return list(reader), fieldNames, smilesColumn


def writeSmilesCsv(outFile, rows, fieldNames):
    '''Write the rows back, keeping the columns of the input csv'''
    with open(outFile, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldNames)
        writer.writeheader()
        writer.writerows(rows)


def checkSubstructures(smilesList, names):
    '''Check that every warhead (all but the last SMILES) is a substructure of the
    PROTAC (the last one), as PROTAC-Model requires'''
    mols = [Chem.MolFromSmiles(smi) for smi in smilesList]

    for smi, mol, name in zip(smilesList, mols, names):
        if mol is None:
            raise ValueError('RDKit cannot parse the normalized SMILES of {}: {}'.
                             format(name, smi))

    protacMol, protacName = mols[-1], names[-1]
    for mol, name in zip(mols[:-1], names[:-1]):
        if not protacMol.HasSubstructMatch(mol):
            raise ValueError(
                'The {} SMILES is not a substructure of the {} SMILES, so PROTAC-Model '
                'cannot match the warhead inside the PROTAC:\n  {}: {}\n  {}: {}'.
                format(name, protacName, name, Chem.MolToSmiles(mol),
                       protacName, Chem.MolToSmiles(protacMol)))


def parseArgs():
    parser = argparse.ArgumentParser(description='Normalize the SMILES of a csv written '
                                                 'by rdkit_IO.py')
    parser.add_argument('-i', '--inputFile', required=True,
                        help='CSV of SMILES to normalize')
    parser.add_argument('-o', '--outputFile', default=None,
                        help='Where to write the normalized csv. Defaults to the input '
                             'file, which is overwritten')
    parser.add_argument('--check-substructures', action='store_true',
                        dest='checkSubstructures',
                        help='Check that every SMILES but the last is a substructure of '
                             'the last one (the PROTAC)')
    parser.add_argument('--quiet', action='store_true',
                        help='Do not print the normalized SMILES')

    return parser.parse_args()


def main():
    args = parseArgs()
    RDLogger.DisableLog('rdApp.*')

    inFile = os.path.abspath(args.inputFile)
    outFile = os.path.abspath(args.outputFile) if args.outputFile else inFile

    if not os.path.isfile(inFile):
        sys.exit('SMILES csv not found: {}'.format(inFile))

    try:
        rows, fieldNames, smilesColumn = readSmilesCsv(inFile)
    except ValueError as e:
        sys.exit(str(e))

    if not rows:
        sys.exit('No SMILES found in {}'.format(inFile))

    names = [row.get('name') or 'molecule {}'.format(i + 1) for i, row in enumerate(rows)]

    try:
        for row, name in zip(rows, names):
            smi = (row[smilesColumn] or '').strip()
            if not smi:
                raise ValueError('Empty SMILES for {} in {}'.format(name, inFile))

            row[smilesColumn] = normalizeSmiles(smi)

        if args.checkSubstructures and len(rows) > 1:
            checkSubstructures([row[smilesColumn] for row in rows], names)

    except (ValueError, Chem.rdchem.MolSanitizeException) as e:
        sys.exit(str(e))

    writeSmilesCsv(outFile, rows, fieldNames)

    if not args.quiet:
        for row, name in zip(rows, names):
            print('{}: {}'.format(name, row[smilesColumn]))
        print('Normalized SMILES csv saved to: {}'.format(outFile))


if __name__ == '__main__':
    main()

