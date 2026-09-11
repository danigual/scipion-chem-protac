#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Python 2 driver for PROTAC-Model (gaoqiweng/PROTAC-Model). Launched as a separate
# process (not imported) by rosetta/protocols/protocol_protac_model.py, with the
# PROTAC_MODEL_PYTHON_HOME conda env activated - see rosetta/__init__.py's
# runCondaScript(). Runs under whatever cwd the caller chose (extra/frodock/ for
# --phase frodock/filter, extra/rosetta/ for --phase refine), matching the layout
# main.py itself builds around fro.frodock()/fro.filter_frodock()/ros.rosetta().
#
# This file only stages input files and calls ONE of PROTAC-Model's own high-level
# functions per phase: splitting those functions further would mean reimplementing their
# internal orchestration (chained os.chdir, shell calls, a multiprocessing.Pool per pose)
# ourselves, which is exactly what calling the original code is meant to avoid.

import argparse
import glob
import os
import shutil
import sys

# PROTAC_MODEL_HOME is set by the caller's extraEnvDict (rosetta/__init__.py's
# PROTAC_MODEL_DIC resolved via getProtacModelScript()) - the repo root containing
# main.py and utils/.
sys.path.insert(0, os.environ['PROTAC_MODEL_HOME'])
import utils.preprocess as pre
import utils.frodock as fro
import utils.rosetta as ros


def _removeIfExists(*paths):
    """ Best-effort delete, used before any phase that (via the original code) writes
    to a file in append mode (>> or 2>>). Without this, retrying a Scipion step would
    double-count lines in these files and corrupt the pose counting/filtering that reads
    them back. """
    for path in paths:
        if os.path.exists(path):
            os.remove(path)


def runFrodock(args):
    """ --phase frodock: replicates what main.py does before calling fro.frodock() -
    chain-ID renaming (pre.alter_pro_chain), staging the PROTAC smiles and the optional
    E3 ligand conformers into cwd - then calls fro.frodock(site). """
    pre.alter_pro_chain(args.receptor, args.target, 'receptor.pdb', 'target.pdb')

    with open('protac.smi', 'w') as f:
        f.write(args.smiles.strip() + '\n')

    if args.e3lig1 and args.e3lig2:
        shutil.copy(args.e3lig1, 'rec_lig_1.sdf')
        shutil.copy(args.e3lig2, 'rec_lig_2.sdf')

    # Fix #4: frodock() ends with 'frodockview ... >> frodock_score.txt' (append).
    _removeIfExists('frodock_score.txt')

    fro.frodock(args.site)


def runFilter(args):
    """ --phase filter: same cwd as --phase frodock (extra/frodock/, already populated
    by it). Calls fro.filter_frodock(), which internally appends to results_voromqa (once
    per pose, from inside a multiprocessing.Pool) and to vina/score_all_top1 /
    vina/score_filter (fix #4 - see _removeIfExists above). """
    _removeIfExists('results_voromqa', 'addH_log')
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    fro.filter_frodock(args.cpu, args.ligLocateNum, args.targetSmi, args.recSmi)


def runRefine(args):
    """ --phase refine: cwd is extra/rosetta/, a sibling of extra/frodock/ (ros.rosetta()
    uses hardcoded relative paths like '../frodock/...'). Stages the files ros.rosetta()
    itself would try to 'cp ... {a,b,c}' via os.system() - that brace expansion is a
    bash-ism and fails when the subprocess runs under dash/sh instead of bash, so it's
    done here with shutil instead - before calling ros.rosetta(). """
    frodockDir = os.path.join('..', 'frodock')
    for src in glob.glob(os.path.join(frodockDir, 'rec_lig_*.sdf')):
        shutil.copy(src, '.')
    shutil.copy(os.path.join(frodockDir, 'rec_lig.sdf'), '.')
    shutil.copy(os.path.join(frodockDir, 'target_lig.sdf'), '.')
    shutil.copy(os.path.join(frodockDir, 'protac.smi'), '.')

    # Fix #4, same reasoning as runFilter above - ros.rosetta()'s own filtering() also
    # appends to results_voromqa/vina score files.
    _removeIfExists('results_voromqa', 'addH_log')
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    ros.rosetta(args.cpu, args.ligLocateNum, args.targetSmi, args.recSmi)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', required=True, choices=['frodock', 'filter', 'refine'])
    # --phase frodock
    parser.add_argument('--receptor')
    parser.add_argument('--target')
    parser.add_argument('--smiles')
    parser.add_argument('--site')
    parser.add_argument('--e3lig1', default=None)
    parser.add_argument('--e3lig2', default=None)
    # --phase filter / refine
    parser.add_argument('--cpu', type=int, default=1)
    parser.add_argument('--lig-locate-num', dest='ligLocateNum', type=int, default=1)
    parser.add_argument('--target-smi', dest='targetSmi', default='none')
    parser.add_argument('--rec-smi', dest='recSmi', default='none')
    return parser.parse_args()


if __name__ == '__main__':
    args = parseArgs()
    if args.phase == 'frodock':
        runFrodock(args)
    elif args.phase == 'filter':
        runFilter(args)
    elif args.phase == 'refine':
        runRefine(args)
