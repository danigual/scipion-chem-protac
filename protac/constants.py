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

# Plugin version. Same split as the rest of the ecosystem: the number lives here and
# __init__.py only does '__version__ = ALPHA_VERSION'.
ALPHA_VERSION = '0.1'

# ------------------------------- Package dictionaries -------------------------------
# FRODOCK is a separate external tool, used by the PROTAC-Model pipeline
# for the initial global protein-protein docking step. 3.12 is the latest stable
# release and the one PROTAC-Model's own README points at.
FRODOCK_DIC = {'name': 'frodock', 'version': '3.12', 'home': 'FRODOCK_HOME'}

# Used later by filterPosesStep (PROTAC-Model's filter_frodock()): ADFRsuite binaries,
# Vina, Voromqa, and FCC's clustering scripts. None of them is license-gated (unlike
# Rosetta below), so these get real defineBinaries() via InstallHelper.
# ADFRsuite's 1.0 is not a free choice: its download URL and the folder name its
# installer extracts (ADFRsuite_x86_64Linux_1.0) both hardcode it.
ADFRSUITE_DIC = {'name': 'adfrsuite', 'version': '1.0', 'home': 'ADFRSUITE_HOME'}
# PROTAC-Model itself pins nothing here: its README still points at the pre-1.2
# vina.scripps.edu download. 1.2.3 is chosen to match scipion-chem-autodock, and the
# runtime patches in run_protac_model.py (uninitialised grid box in --score_only,
# reworded score line) were checked against Vina's own source at 1.2.2, 1.2.3 and 1.2.5.
VINA_DIC = {'name': 'vina', 'version': '1.2.3', 'home': 'VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'version': '1.29.4816', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'version': 'latest', 'home': 'FCC_HOME'}

# PROTAC-Model's own code (main.py, utils/*) is called directly, not reimplemented
# Its code is genuine Python 2, so PROTAC_MODEL_PYTHON_HOME is a dedicated
# Python 2.7+rdkit conda env, never scipion3's own.
PROTAC_MODEL_DIC = {'name': 'protac-model', 'version': 'latest', 'home': 'PROTAC_MODEL_HOME'}
PROTAC_MODEL_PYTHON_DIC = {'name': 'protac-model-python', 'version': '2.7',
                           'home': 'PROTAC_MODEL_PYTHON_HOME'}

# The 4 programs PROTAC-Model's own utils/frodock.py resolves, each shipped as an
# intel/gcc build pair ('<name>' / '<name>_gcc').
FRODOCK_BINARIES = ['frodockgrid', 'frodock', 'frodockcluster', 'frodockview']

# ---------------------------- ProtPROTACTransplantWarhead ----------------------------
# A: a source residue counts as "pocket" if it has an atom within this distance of the
# warhead being transplanted.
POCKET_CUTOFF = 6.0
# A: ligand-protein contacts closer than this, in the target, are counted as clashes.
CLASH_CUTOFF = 2.0
# Minimum number of source pocket residues that must have an aligned equivalent in the
# target for the pocket-restricted superposition to be considered reliable.
MIN_POCKET_PAIRS = 10
# Canonical (uppercase) one-letter code for residues that appear as HETATM in a PDB but
# are actually part of the main chain. Deliberately NOT pwchem.utils.MODIFIED_RESIDUES3TO1:
# that dict uses arbitrary lowercase codes for a different purpose (mutation wizards) that
# fall outside BLOSUM62's alphabet and would silently break sequence alignment here.
MODRES_TO_CANONICAL = {'PTR': 'Y', 'TPO': 'T', 'SEP': 'S', 'MSE': 'M'}
