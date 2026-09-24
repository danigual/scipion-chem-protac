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

ALPHA_VERSION = '0.1'

# ------------------------------- Package dictionaries -------------------------------
# Global protein-protein docking for PROTAC-Model. 3.12 is the release its README uses.
FRODOCK_DIC = {'name': 'frodock', 'version': '3.12', 'home': 'FRODOCK_HOME'}

# Tools used by PROTAC-Model's pose filtering. ADFRsuite 1.0 is fixed by its download
# URL and installer folder name.
ADFRSUITE_DIC = {'name': 'adfrsuite', 'version': '1.0', 'home': 'ADFRSUITE_HOME'}
# Same version as scipion-chem-autodock. The Vina fixes in run_protac_model.py were
# checked against Vina's source for 1.2.2, 1.2.3 and 1.2.5. Not VINA_HOME: autodock uses
# that name for a tree without bin/vina.
VINA_DIC = {'name': 'vina', 'version': '1.2.3', 'home': 'PROTAC_VINA_HOME'}
VOROMQA_DIC = {'name': 'voromqa', 'version': '1.29.4816', 'home': 'VOROMQA_HOME'}
FCC_DIC = {'name': 'fcc', 'version': 'latest', 'home': 'FCC_HOME'}

# PROTAC-Model is Python 2, so it gets its own Python 2.7 + RDKit env.
PROTAC_MODEL_DIC = {'name': 'protac-model', 'version': 'latest', 'home': 'PROTAC_MODEL_HOME'}
PROTAC_MODEL_PYTHON_DIC = {'name': 'protac-model-python', 'version': '2.7',
                           'home': 'PROTAC_MODEL_PYTHON_HOME'}

# The FRODOCK programs PROTAC-Model calls, each shipped as '<name>' (intel) and
# '<name>_gcc'.
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
# One-letter codes for modified residues often stored as HETATM. pwchem's
# MODIFIED_RESIDUES3TO1 uses lowercase codes that BLOSUM62 can't align.
MODRES_TO_CANONICAL = {'PTR': 'Y', 'TPO': 'T', 'SEP': 'S', 'MSE': 'M'}

# ---------------------------------- ProtPRosettaC ----------------------------------
PROSETTAC_DIC = {'name': 'prosettac', 'version': 'latest', 'home': 'PROSETTAC_HOME'}

# XXX: rdkit/numpy/scikit-learn versions not pinned yet.
# TODO: set the installed versions (conda-forge). The placeholder ends up in the env
# name, so an env installed now will need reinstalling under the final one.
PROSETTAC_PYTHON_DIC = {'name': 'prosettac-python', 'version': 'XXX',
                        'home': 'PROSETTAC_PYTHON_HOME'}

# Bare Python 2.7 for Rosetta's molfile_to_params.py, which only uses the standard
# library. Kept apart from PROTAC-Model's env to avoid its 2016 RDKit pin.
PROSETTAC_PYTHON2_DIC = {'name': 'prosettac-python2', 'version': '2.7',
                         'home': 'PROSETTAC_PYTHON2_HOME'}

# Installed by hand (academic license), like Rosetta.
# TODO: set the installed version (https://bioinfo3d.cs.tau.ac.il/PatchDock/).
PATCHDOCK_DIC = {'name': 'patchdock', 'version': 'XXX', 'home': 'PATCHDOCK_HOME'}

# rosetta_scripts builds, in the order the shim looks for them. PRosettaC expects
# .default. (a source build); prebuilt bundles ship .static., MPI builds .mpi.
ROSETTA_SCRIPTS_BUILDS = ['default', 'static', 'mpi']
