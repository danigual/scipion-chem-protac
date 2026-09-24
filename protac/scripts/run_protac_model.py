#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Python 2 driver for PROTAC-Model, run by ProtPROTACModel inside the
# PROTAC_MODEL_PYTHON_HOME env, from extra/frodock/ (frodock, filter) or extra/rosetta/
# (refine). Each phase stages its inputs and calls one of PROTAC-Model's own
# high-level functions; splitting them further would mean reimplementing their
# internal orchestration.

import argparse
import functools
import glob
import math
import os
import re
import shutil
import sys
import traceback

# Repo root of the PROTAC-Model checkout, set by the protocol.
sys.path.insert(0, os.environ['PROTAC_MODEL_HOME'])
import utils.preprocess as pre
import utils.frodock as fro
import utils.rosetta as ros


# ----------------------------- Robustness shims ------------------------------
# PROTAC-Model assumes every shell pipeline it launches prints something and that no pose
# ever fails. preprocess.obenergy_vina() reads each conformer's Vina score with
# os.popen(...).read().splitlines()[0], so one conformer without output raises IndexError
# inside a Pool worker and kills the whole step. Rather than patching the checkout on
# disk, the driver patches the loaded modules at runtime. The patches reach the Pool
# workers because they are applied before the Pool is forked, and both patched names
# (the module's 'os' and 'filtering') are looked up at call time.

_MAX_WARNINGS_PER_KIND = 20
_warningCounts = {}

# One line per pose skipped by the fail-safe filtering(), read by _reportPoseFailures().
POSE_FAILURES_LOG = 'protac_pose_failures.log'


def _warn(message):
    """ To stdout, next to the original code's output. One flushed write, so messages
    from different Pool workers don't interleave. """
    sys.stdout.write(message)
    sys.stdout.flush()


def _warnCapped(kind, message):
    """ Capped per kind and process: a systematic problem would otherwise warn once per
    conformer of every pose. """
    count = _warningCounts.get(kind, 0) + 1
    _warningCounts[kind] = count
    if count > _MAX_WARNINGS_PER_KIND:
        return
    if count == _MAX_WARNINGS_PER_KIND:
        message += ('[protac] WARNING: further "%s" warnings from this process are '
                    'suppressed.\n' % kind)
    _warn(message)


class _GuardedPipe(object):
    """ Stand-in for os.popen()'s file object, holding output already read. PROTAC-Model
    only uses read() and 'with' on it. """

    def __init__(self, content):
        self._content = content

    def read(self):
        return self._content

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, excTb):
        return False


# --------------------- Vina --score_only: grid box + score parsing ---------------------
# obenergy_vina() scores each conformer with
#     vina --score_only --receptor ... --ligand ... | grep Affinity | cut -d" " -f2
# which was written for pre-1.2 Vina and breaks on current versions in two ways:
#   1. It passes no search space. Vina 1.2 then builds its maps from uninitialised
#      values, so every conformer gets a nonsense energy or "ligand is outside the grid
#      box", and the step later dies on a missing vina/score_all_top1.
#   2. Vina 1.2 no longer prints an "Affinity:" line; the score line is now
#      "Estimated Free Energy of Binding : ...".
# The popen() shim recognises this command, reruns Vina with a box around the ligand and
# returns the single score the grep|cut tail was meant to produce.

# Vina's own --autobox padding. Kept minimal: map cost grows with box volume, and this
# runs for every conformer of every pose.
_VINA_BOX_PADDING = 4.0

# Returned instead of a score that could not be obtained (see _GuardedOs).
_ZERO_SCORE = '0\n'

# '\S+' stops before the '|' of the grep pipeline, so it captures just the path.
_VINA_LIGAND_RE = re.compile(r'--ligand(?:\s+|=)(\S+)')

# Both wordings of the score line (Vina 1.1.x and 1.2.x), anchored at line start so the
# per-term breakdown lines of 1.2.x never match.
_VINA_SCORE_RE = re.compile(
    r'^[ \t]*(?:Affinity|Estimated Free Energy of Binding)[ \t]*:[ \t]*'
    r'([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)',
    re.MULTILINE)


def _readPdbqtCoords(pdbqtPath):
    """ Every (x, y, z) of the ATOM/HETATM records of a PDBQT file (standard PDB
    coordinate columns). Malformed lines are skipped rather than read as 0.0, which would
    stretch the box towards the origin. """
    coords = []
    with open(pdbqtPath) as pdbqtFile:
        for line in pdbqtFile:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return coords


def _vinaBoxFlags(ligandPath):
    """ Vina box flags for a box around the ligand with _VINA_BOX_PADDING on every face,
    or None if the file has no usable atom. Same size formula as Vina's autobox, but
    centred on the bounding box rather than the centroid, which gives a smaller box for
    elongated molecules like a PROTAC. """
    try:
        coords = _readPdbqtCoords(ligandPath)
    except (IOError, OSError):
        return None
    if not coords:
        return None

    centers, sizes = [], []
    for axis in range(3):
        values = [coord[axis] for coord in coords]
        low, high = min(values), max(values)
        centers.append((low + high) / 2.0)
        sizes.append(math.ceil((high - low) / 2.0 + _VINA_BOX_PADDING) * 2.0)
    return ('--center_x %.3f --center_y %.3f --center_z %.3f '
            '--size_x %.3f --size_y %.3f --size_z %.3f'
            % (centers[0], centers[1], centers[2], sizes[0], sizes[1], sizes[2]))


def _isVinaScoreOnly(cmd):
    """ True only for obenergy_vina()'s Vina call, and only while it still lacks a
    search space of its own. """
    if '--score_only' not in cmd:
        return False
    return not ('--center_x' in cmd or '--maps' in cmd or '--autobox' in cmd)


def _runVinaScoreOnly(realOs, cmd, args, kwargs):
    """ Reruns the Vina part of cmd (before the first '|') with a box and returns one
    line with the score, or '0' on failure. None if cmd isn't shaped as expected. """
    ligandMatch = _VINA_LIGAND_RE.search(cmd)
    if ligandMatch is None or not ligandMatch.group(1).endswith('.pdbqt'):
        return None
    ligandPath = ligandMatch.group(1)

    boxFlags = _vinaBoxFlags(ligandPath)
    if boxFlags is None:
        _warnCapped('vina-no-box',
                    '[protac] WARNING: cannot read ligand coordinates from %s, skipping '
                    'this conformer (scoring it as 0).\n' % ligandPath)
        return _ZERO_SCORE

    # 2>&1 so a Vina error ends up in the warning below.
    vinaCmd = '%s %s 2>&1' % (cmd.split('|', 1)[0].strip(), boxFlags)
    pipe = realOs.popen(vinaCmd, *args, **kwargs)
    try:
        output = pipe.read()
    finally:
        status = pipe.close()
    if status and status >> 8 in (126, 127):
        raise _FatalError('Vina cannot be run: %s' % output.strip())

    scoreMatch = _VINA_SCORE_RE.search(output)
    if scoreMatch is None:
        _warnCapped('vina-no-score',
                    '[protac] WARNING: no score in the output of "%s", scoring this '
                    'conformer as 0. Vina said:\n%s\n' % (vinaCmd, output.strip()))
        return _ZERO_SCORE
    return '%s\n' % scoreMatch.group(1)


class _GuardedOs(object):
    """ The 'os' module as seen from utils/preprocess.py, except that popen() fixes the
    Vina score call and never returns empty output. '0' is a safe stand-in for both of
    preprocess.py's popen() calls: a 0 score drops that conformer (awk '$2<0') without
    shifting the one-line-per-conformer files it is pasted against, and a 0 count drops
    the pose. """

    def __init__(self, realOs):
        self._os = realOs

    def __getattr__(self, name):
        return getattr(self._os, name)

    def popen(self, cmd, *args, **kwargs):
        if _isVinaScoreOnly(cmd):
            # None: not the expected shape after all, use the generic path.
            content = _runVinaScoreOnly(self._os, cmd, args, kwargs)
            if content is not None:
                return _GuardedPipe(content)

        pipe = self._os.popen(cmd, *args, **kwargs)
        try:
            content = pipe.read()
        finally:
            pipe.close()
        lines = content.splitlines()
        if not lines or not lines[0].strip():
            _warnCapped('empty-output',
                        '[protac] WARNING: no output from shell command, using 0 instead: '
                        '%s\n' % cmd)
            content = _ZERO_SCORE
        return _GuardedPipe(content)


class _FatalError(Exception):
    """ Not a per-pose failure: stops the whole run. """


def _makeFailSafeFiltering(module, poseIndex):
    """ Wraps utils.<module>.filtering() so a failing pose is logged and skipped instead
    of killing the Pool. Safe because a pose only reaches results_voromqa as its last
    action, and everything after the Pool reads from there. functools.wraps is required:
    pool.map() pickles the callable by its qualified name. poseIndex is where the pose id
    sits in para_list (0 in utils/frodock.py, 4 in utils/rosetta.py). """
    original = module.filtering

    @functools.wraps(original)
    def filtering(paraList):
        try:
            return original(paraList)
        except _FatalError:
            raise
        except Exception:
            try:
                poseId = str(paraList[poseIndex])
            except Exception:
                poseId = '<unknown>'
            _warn('[protac] WARNING: skipping pose %s, it failed inside PROTAC-Model\'s '
                  'filtering():\n%s' % (poseId, traceback.format_exc()))
            try:
                with open(POSE_FAILURES_LOG, 'a') as logFile:
                    logFile.write('%s\n' % poseId)
            except Exception:
                pass
        return None

    return filtering


def _forceCLocale():
    """ The original code filters numbers with awk. Under a locale with a decimal comma
    (es_ES, for instance) mawk compares '120.5' as text, every conformer is dropped and
    the step dies later with no clear cause. Child processes inherit this C locale. """
    os.environ['LC_ALL'] = 'C'


def installRuntimeFixes():
    """ Called before any phase, so the patches are in place before a Pool is forked. """
    _forceCLocale()
    if not isinstance(pre.os, _GuardedOs):
        pre.os = _GuardedOs(pre.os)
    fro.filtering = _makeFailSafeFiltering(fro, 0)
    ros.filtering = _makeFailSafeFiltering(ros, 4)


def _reportPoseFailures():
    """ Summary of the skipped poses. Called from a finally block: if every pose failed,
    the original code dies later on a missing file, and this line explains why. """
    if not os.path.exists(POSE_FAILURES_LOG):
        return
    with open(POSE_FAILURES_LOG) as logFile:
        failures = [line.strip() for line in logFile.read().splitlines() if line.strip()]
    if not failures:
        return
    shown = failures[:50]
    if len(failures) > len(shown):
        shown.append('...')
    _warn('[protac] WARNING: %d pose(s) skipped because they failed inside PROTAC-Model: '
          '%s\n' % (len(failures), ', '.join(shown)))


def _removeIfExists(*paths):
    """ Clears files the original code appends to, so a retried step doesn't count the
    same lines twice. """
    for path in paths:
        if os.path.exists(path):
            os.remove(path)


def runFrodock(args):
    """ Does what main.py does before fro.frodock(): renames chains, stages the PROTAC
    SMILES and the optional E3 ligand conformers. """
    pre.alter_pro_chain(args.receptor, args.target, 'receptor.pdb', 'target.pdb')

    with open('protac.smi', 'w') as f:
        f.write(args.smiles.strip() + '\n')

    if args.e3lig1 and args.e3lig2:
        shutil.copy(args.e3lig1, 'rec_lig_1.sdf')
        shutil.copy(args.e3lig2, 'rec_lig_2.sdf')

    # frodock is called without --soap, so it expects soap.bin in cwd.
    shutil.copy(os.path.join(os.environ['FRODOCK'], 'bin', 'soap.bin'), 'soap.bin')

    # frodock() appends to it.
    _removeIfExists('frodock_score.txt')

    fro.frodock(args.site)


def runFilter(args):
    """ Runs in the frodock phase's directory. filter_frodock() appends to
    results_voromqa and the vina score files, so those are cleared first. """
    _removeIfExists('results_voromqa', 'addH_log', POSE_FAILURES_LOG)
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    try:
        fro.filter_frodock(args.cpu, args.ligLocateNum, args.targetSmi, args.recSmi)
    finally:
        _reportPoseFailures()


def runRefine(args):
    """ Runs in extra/rosetta/, next to extra/frodock/ (ros.rosetta() uses '../frodock/'
    paths). Copies the files ros.rosetta() would copy with a bash brace expansion, which
    fails under sh. """
    frodockDir = os.path.join('..', 'frodock')
    for src in glob.glob(os.path.join(frodockDir, 'rec_lig_*.sdf')):
        shutil.copy(src, '.')
    shutil.copy(os.path.join(frodockDir, 'rec_lig.sdf'), '.')
    shutil.copy(os.path.join(frodockDir, 'target_lig.sdf'), '.')
    shutil.copy(os.path.join(frodockDir, 'protac.smi'), '.')

    _removeIfExists('results_voromqa', 'addH_log', POSE_FAILURES_LOG)
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    try:
        ros.rosetta(args.cpu, args.ligLocateNum, args.targetSmi, args.recSmi)
    finally:
        _reportPoseFailures()


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', required=True, choices=['frodock', 'filter', 'refine'])
    parser.add_argument('--receptor')
    parser.add_argument('--target')
    parser.add_argument('--smiles')
    parser.add_argument('--site')
    parser.add_argument('--e3lig1', default=None)
    parser.add_argument('--e3lig2', default=None)
    parser.add_argument('--cpu', type=int, default=1)
    parser.add_argument('--lig-locate-num', dest='ligLocateNum', type=int, default=1)
    parser.add_argument('--target-smi', dest='targetSmi', default='none')
    parser.add_argument('--rec-smi', dest='recSmi', default='none')
    return parser.parse_args()


if __name__ == '__main__':
    args = parseArgs()
    installRuntimeFixes()
    if args.phase == 'frodock':
        runFrodock(args)
    elif args.phase == 'filter':
        runFilter(args)
    elif args.phase == 'refine':
        runRefine(args)
