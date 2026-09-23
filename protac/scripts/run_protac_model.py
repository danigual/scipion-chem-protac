#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Python 2 driver for PROTAC-Model (gaoqiweng/PROTAC-Model). Launched as a separate
# process by protac/protocols/protocol_protac_model.py under the PROTAC_MODEL_PYTHON_HOME
# conda env (see protac/__init__.py's runCondaScript()), from whatever cwd the caller
# chose (extra/frodock/ for --phase frodock/filter, extra/rosetta/ for --phase refine).
#
# Only stages inputs and calls ONE of PROTAC-Model's own high-level functions per phase;
# splitting them further would mean reimplementing their internal orchestration.

import argparse
import functools
import glob
import math
import os
import re
import shutil
import sys
import traceback

# PROTAC_MODEL_HOME is set by the caller's extraEnvDict (rosetta/__init__.py's
# PROTAC_MODEL_DIC resolved via getProtacModelScript()) - the repo root containing
# main.py and utils/.
sys.path.insert(0, os.environ['PROTAC_MODEL_HOME'])
import utils.preprocess as pre
import utils.frodock as fro
import utils.rosetta as ros


# ----------------------------- Robustness shims ------------------------------
# PROTAC-Model's own code assumes every shell pipeline it launches produces output, and
# that no pose ever fails. The worst offender is preprocess.obenergy_vina(), which does
#
#     with os.popen('... vina --score_only ... | grep Affinity | cut -d" " -f2') as f:
#         score = f.read().splitlines()[0]
#
# once per conformer (up to 100) of every pose (hundreds of them). A single conformer for
# which vina prints no "Affinity" line - for any reason: a transient prepare_ligand
# failure, resource contention between Pool workers, a weird conformer geometry... - makes
# that [0] raise IndexError inside a multiprocessing.Pool worker; pool.map() re-raises it
# in the parent and the whole step dies, throwing away hours of docking because of one
# conformer. (That command is also wrong in two further, systematic ways against any Vina
# we can install today - no grid box, and the wrong wording to grep for. Both are handled
# in the "Vina --score_only" section below, which is why the guard here is now only the
# last-resort net it was meant to be.)
#
# We neither reimplement that logic nor patch
# the PROTAC-Model checkout on disk (PROTAC_MODEL_HOME may point at a clone this plugin
# never made, and a reinstall would silently drop the patch). Instead we harden it at
# runtime from here, and both patches survive into the Pool workers because:
#   - Python 2's multiprocessing forks on Linux, so a worker inherits the parent's already
#     patched module objects, as long as the patch is applied BEFORE the Pool is created
#     (i.e. before calling fro.filter_frodock() / ros.rosetta(), which build their own);
#   - both patched names are resolved dynamically at call time, not captured beforehand:
#     obenergy_vina() reaches popen through its module global 'os', and pool.map() pickles
#     'filtering' by qualified name (utils.frodock.filtering), which the worker resolves
#     against its inherited copy of the module.

_MAX_WARNINGS_PER_KIND = 20
_warningCounts = {}

# Appended to (one line per failed pose) by the fail-safe filtering() below, read back by
# _reportPoseFailures(). Lives in the phase's working directory, like everything else the
# original code writes.
POSE_FAILURES_LOG = 'protac_pose_failures.log'


def _warn(message):
    """ Warnings go to stdout (not stderr) so they show up in Scipion's run.stdout next to
    the output of the original code, in order. Written in one call and flushed: Pool
    workers share this fd, and a single short write to a pipe is not interleaved. """
    sys.stdout.write(message)
    sys.stdout.flush()


def _warnCapped(kind, message):
    """ Capped per kind and per process (each forked worker gets its own counters): when
    something systematic is broken these fire once per conformer - up to 100 conformers x
    hundreds of poses - and an uncapped log would be tens of MB. """
    count = _warningCounts.get(kind, 0) + 1
    _warningCounts[kind] = count
    if count > _MAX_WARNINGS_PER_KIND:
        return
    if count == _MAX_WARNINGS_PER_KIND:
        message += ('[protac] WARNING: further "%s" warnings from this process are '
                    'suppressed.\n' % kind)
    _warn(message)


class _GuardedPipe(object):
    """ Stand-in for the file object os.popen() returns, holding output already read. Only
    implements what PROTAC-Model actually uses on it (read(), and the 'with' protocol). """

    def __init__(self, content):
        self._content = content

    def read(self):
        return self._content

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, excTb):
        return False


# --------------------- Vina --score_only: grid box + score parsing ---------------------
# PROTAC-Model was written against pre-1.2 AutoDock Vina (its README points at the old
# vina.scripps.edu download), and obenergy_vina() scores each conformer with
#
#     $VINA/bin/vina --score_only --receptor <...>.pdbqt --ligand <...>.pdbqt \
#         | grep Affinity | cut -d" " -f2
#
# Two things about that command are broken against every Vina we can actually install
# today (both verified against ccsb-scripps/AutoDock-Vina, tags v1.2.2 and v1.2.5):
#
#   1. No search space. In src/main/main.cpp the score_only branch is
#          if ((score_only || local_only) && autobox) { ...from ligand... }
#          else v.compute_vina_maps(center_x, ..., size_z, grid_spacing, force_even_voxels);
#      and center_x/size_x are plain 'double center_x;' with no initialiser and no
#      vm.count() check. So without a box Vina silently builds its affinity maps from
#      uninitialised stack garbage (e.g. "Center: X 1.58101e-322 ...")
#      - and either returns a nonsense energy (+2.3e9 kcal/mol, which the
#      downstream awk '$2<0' filter then discards) or dies in Vina::score() on
#      m_grid.is_in_grid(m_model) with "The ligand is outside the grid box". Either way no
#      conformer of any pose ever survives, vina/score_all_top1 is never written, and
#      filter_frodock() dies much later on IOError: 'vina/score_all_top1'. This is not a
#      version regression we can pin our way out of (1.2.2 was already tried): the box is
#      simply required, and this build has no --autobox either.
#   2. No "Affinity" line to grep. Vina::show_score() in src/lib/vina.cpp prints
#      "Estimated Free Energy of Binding   : <x> (kcal/mol) [=(1)+(2)+(3)+(4)]"; the
#      "Affinity: <x> (kcal/mol)" line the pipeline greps for is the pre-1.2 wording. So
#      even with a correct box the grep|cut pipeline would yield nothing.
#
# Both are fixed here, inside the popen() shim we already own, by taking over the whole
# command when it is recognisably this one: we compute a real box from the ligand PDBQT,
# re-run Vina with it, and hand back the single number the original grep|cut pipeline was
# meant to produce. No PROTAC-Model logic is reimplemented - obenergy_vina() still decides
# what to do with the score.

# Padding added around the ligand on every face, in Angstrom. 4.0 is Vina's own --autobox
# buffer ('double buffer_size = 4;' in src/main/main.cpp), i.e. what Vina considers enough
# clearance to score a ligand. Kept at exactly that and no more on purpose: Vina builds
# real affinity maps for --score_only, so map construction cost grows with the box volume
# and this runs once per conformer (up to 100) of every pose.
_VINA_BOX_PADDING = 4.0

# What we hand back in place of a score we could not obtain. See _GuardedOs' docstring for
# why 0 (rather than no line at all) is the safe stand-in here.
_ZERO_SCORE = '0\n'

# The --ligand argument of the command being intercepted. '\S+' stops at the space before
# the '|' of the grep pipeline, so it captures just the path.
_VINA_LIGAND_RE = re.compile(r'--ligand(?:\s+|=)(\S+)')

# Accepts both wordings of the score line so this keeps working if PROTAC_MODEL's original
# pre-1.2 Vina is ever the one installed:
#   "Affinity: -8.51139 (kcal/mol)"                                     (Vina 1.1.x)
#   "Estimated Free Energy of Binding   : -8.512 (kcal/mol) [=(1)+...]" (Vina 1.2.x)
# Anchored at the start of a line so the "(1) Final Intermolecular Energy : ..." breakdown
# lines 1.2.x prints right after it can never match.
_VINA_SCORE_RE = re.compile(
    r'^[ \t]*(?:Affinity|Estimated Free Energy of Binding)[ \t]*:[ \t]*'
    r'([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)',
    re.MULTILINE)


def _readPdbqtCoords(pdbqtPath):
    """ Every (x, y, z) in a PDBQT file, as float triples.

    Column ranges are not guessed: they are the ones Vina's own reader uses in
    src/lib/parse_pdbqt.cpp,

        vec coords(checked_convert_substring<fl>(str, 31, 38, "Coordinate"),
                   checked_convert_substring<fl>(str, 39, 46, "Coordinate"),
                   checked_convert_substring<fl>(str, 47, 54, "Coordinate"));

    which are 1-indexed and inclusive, hence the [30:38] / [38:46] / [46:54] slices below -
    the same ones PROTAC-Model itself uses on PDB files (see preprocess.lig_around_residue).
    PDBQT only adds the partial charge (cols 69-76) and AutoDock atom type (cols 78-79)
    after the B-factor, so the coordinate columns are the plain PDB ones.

    Everything that is not an ATOM/HETATM record (REMARK, ROOT/BRANCH/TORSDOF, MODEL...)
    carries no coordinates and is skipped, as are records too short or malformed to parse -
    a truncated line must not silently contribute a 0.0 and blow the box up towards the
    origin. """
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
    """ '--center_x ... --size_z ...' for a box that contains ligandPath with
    _VINA_BOX_PADDING of clearance on every face, or None if the file cannot be read or
    holds no usable atom.

    The size formula mirrors Vina::grid_dimensions_from_ligand() (src/lib/vina.cpp),
    'std::ceil((max_distance[j] + buffer_size) * 2)', with the box centred on the ligand's
    bounding box instead of on its centroid: that is the same guarantee (>= padding of
    clearance on every face) with a strictly smaller box, since the centroid can sit much
    closer to one face than the other for an elongated molecule - and a PROTAC is about as
    elongated as small molecules get. """
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
    """ True only for obenergy_vina()'s Vina call. preprocess.py's other popen() - an
    'awk ... | wc -l' over obenergy_filter_<pose> - has no --score_only and must keep going
    through the generic path untouched.

    A command that already carries a search space (--center_x), precomputed maps (--maps)
    or --autobox is left alone too, so this stays a no-op if PROTAC-Model ever fixes its
    own call. """
    if '--score_only' not in cmd:
        return False
    return not ('--center_x' in cmd or '--maps' in cmd or '--autobox' in cmd)


def _runVinaScoreOnly(realOs, cmd, args, kwargs):
    """ Runs obenergy_vina()'s Vina call with a real grid box and returns exactly what its
    'grep Affinity | cut -d" " -f2' pipeline was supposed to return: one line holding the
    score. Never returns empty output - '0' on any failure, which obenergy_vina() handles
    as "this conformer scored nothing" (see _GuardedOs' docstring).

    Only the part of cmd before the first '|' is re-run: that is the Vina invocation, and
    the grep|cut tail is replaced by _VINA_SCORE_RE. Splitting on '|' is safe for this
    specific command - the only quoting in it is cut's -d" ", which contains no pipe. """
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

    # 2>&1 so a Vina error message ends up in the text we report below instead of being
    # scattered across Scipion's log out of order.
    vinaCmd = '%s %s 2>&1' % (cmd.split('|', 1)[0].strip(), boxFlags)
    pipe = realOs.popen(vinaCmd, *args, **kwargs)
    try:
        output = pipe.read()
    finally:
        pipe.close()

    scoreMatch = _VINA_SCORE_RE.search(output)
    if scoreMatch is None:
        _warnCapped('vina-no-score',
                    '[protac] WARNING: no score in the output of "%s", scoring this '
                    'conformer as 0. Vina said:\n%s\n' % (vinaCmd, output.strip()))
        return _ZERO_SCORE
    return '%s\n' % scoreMatch.group(1)


class _GuardedOs(object):
    """ Drop-in replacement for the 'os' module as seen from utils/preprocess.py. Every
    attribute is the real one except popen(), which (a) gives obenergy_vina()'s Vina call
    the search space it never passes and reads the score back itself (see the block comment
    above), and (b) never returns empty output.

    Why '0' is the right stand-in value for both of preprocess.py's popen call sites:
      - the vina score (its line 37): obenergy_vina writes it to score_<pose> and later
        filters those lines with awk '$2<0', so a 0 drops that conformer. Note we return
        '0' rather than skipping the line altogether: score_<pose> is pasted column-wise
        against obenergy_process_<pose>, which has exactly one line per conformer, so
        dropping a line would shift every later conformer's energy onto the wrong score.
      - the conformer count (its line 49): read back as int(num) > 0, so a 0 drops the
        whole pose.
    In both cases the outcome is "as if this conformer/pose had scored nothing", a state
    the original code already handles (it is what happens when conetnt_score stays empty),
    instead of an IndexError that kills every pose being processed in parallel. """

    def __init__(self, realOs):
        self._os = realOs

    def __getattr__(self, name):
        return getattr(self._os, name)

    def popen(self, cmd, *args, **kwargs):
        if _isVinaScoreOnly(cmd):
            # None means the command was not shaped the way we expect after all; fall
            # through to the generic path rather than guessing.
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


def _makeFailSafeFiltering(module, poseIndex):
    """ Wraps utils.<module>.filtering() so that one failing pose is logged and skipped
    instead of taking down the whole Pool (see the block comment above). Skipping is safe
    by construction: what a pose contributes is appended to results_voromqa as its very
    last action, and the code that runs after the Pool only iterates over the poses that
    made it into results_voromqa - a skipped pose is indistinguishable from one that was
    filtered out on its merits.

    functools.wraps is load-bearing here, not cosmetic: pool.map() pickles the callable by
    qualified name and refuses to if getattr(<its __module__>, <its __name__>) is not the
    object itself, so the wrapper has to claim the very name we install it under.

    poseIndex is where the pose id sits in filtering()'s para_list, which differs between
    utils/frodock.py (0, the pose number) and utils/rosetta.py (4, the model pdb name). """
    original = module.filtering

    @functools.wraps(original)
    def filtering(paraList):
        try:
            return original(paraList)
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
    """ Everything the original code shells out to (awk, sort, wc, obabel, vina...)
    inherits this process' environment, and a lot of it parses numbers with awk. Under a
    locale whose decimal separator is not '.' - es_ES.UTF-8, say - mawk does not recognise
    '120.5' as a number, so obenergy_vina's filter

        awk '{if($2<0 && $4<10000) print $1" "$2" "$4}' obenergy_merge_<pose>

    silently degrades into a *string* comparison for $4, and "120.5" < "10000" is false as
    text. Every conformer of every pose gets dropped, nothing reaches voromqa, and
    filter_frodock dies much later on a missing vina/score_all_top1 with no hint as to why.
    Verified here against mawk 1.3.4 with LC_NUMERIC=es_ES.UTF-8 (0 lines kept) vs LC_ALL=C
    (the expected lines kept). The C locale is what PROTAC-Model implicitly assumes, so we
    pin it instead of depending on how the machine running Scipion happens to be set up.
    Python itself never calls setlocale(), so this only affects the child processes. """
    os.environ['LC_ALL'] = 'C'


def installRuntimeFixes():
    """ Called once before dispatching to a phase, so everything is in place before any
    multiprocessing.Pool is forked. """
    _forceCLocale()
    if not isinstance(pre.os, _GuardedOs):
        pre.os = _GuardedOs(pre.os)
    fro.filtering = _makeFailSafeFiltering(fro, 0)
    ros.filtering = _makeFailSafeFiltering(ros, 4)


def _reportPoseFailures():
    """ Summary of what the fail-safe filtering() skipped. Called from a finally block so
    it prints even when the original code goes on to fail downstream: if every pose was
    skipped, filter_frodock() dies later on a missing vina/score_all_top1, and this is the
    line that explains why. """
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

    # fro.frodock() shells out to 'frodock' without -s/--soap, so it looks for
    # soap.bin in cwd - FRODOCK ships it under bin/, not cwd, and errors out if missing.
    shutil.copy(os.path.join(os.environ['FRODOCK'], 'bin', 'soap.bin'), 'soap.bin')

    # frodock() ends with 'frodockview ... >> frodock_score.txt' (append).
    _removeIfExists('frodock_score.txt')

    fro.frodock(args.site)


def runFilter(args):
    """ --phase filter: same cwd as --phase frodock (extra/frodock/, already populated
    by it). Calls fro.filter_frodock(), which internally appends to results_voromqa (once
    per pose, from inside a multiprocessing.Pool) and to vina/score_all_top1 /
    vina/score_filter. """
    _removeIfExists('results_voromqa', 'addH_log', POSE_FAILURES_LOG)
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    # See "Robustness shims" above - poses that fail are now skipped and reported
    # instead of killing the step from inside a Pool worker.
    try:
        fro.filter_frodock(args.cpu, args.ligLocateNum, args.targetSmi, args.recSmi)
    finally:
        _reportPoseFailures()


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

    # Same reasoning as runFilter above - ros.rosetta()'s own filtering() also
    # appends to results_voromqa/vina score files.
    _removeIfExists('results_voromqa', 'addH_log', POSE_FAILURES_LOG)
    if args.ligLocateNum > 1:
        _removeIfExists(os.path.join('rec_lig_1', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_1', 'vina', 'score_filter'),
                        os.path.join('rec_lig_2', 'vina', 'score_all_top1'),
                        os.path.join('rec_lig_2', 'vina', 'score_filter'))
    else:
        _removeIfExists(os.path.join('vina', 'score_all_top1'),
                        os.path.join('vina', 'score_filter'))

    # Same as runFilter - ros.rosetta() runs the same pose-level Pool.
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
