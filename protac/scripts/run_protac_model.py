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
import os
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
# conformer. Pinning the Vina version
# removed one systematic cause of that, not the fragility of the pattern itself.
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

_MAX_EMPTY_OUTPUT_WARNINGS = 20
_emptyOutputWarnings = 0

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


class _GuardedPipe(object):
    """ Stand-in for the file object os.popen() returns, holding output already read. Only
    implements what PROTAC-Model actually uses on it (read(), and the 'with' protocol). """

    def __init__(self, content):
        self._content = content

    def read(self):
        return self._content

    def readlines(self):
        return self._content.splitlines(True)

    def __iter__(self):
        return iter(self.readlines())

    def close(self):
        return None

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, excTb):
        return False


class _GuardedOs(object):
    """ Drop-in replacement for the 'os' module as seen from utils/preprocess.py. Every
    attribute is the real one except popen(), which never returns empty output.

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

    _FALLBACK = '0\n'

    def __init__(self, realOs):
        self._os = realOs

    def __getattr__(self, name):
        return getattr(self._os, name)

    def popen(self, cmd, *args, **kwargs):
        pipe = self._os.popen(cmd, *args, **kwargs)
        try:
            content = pipe.read()
        finally:
            pipe.close()
        lines = content.splitlines()
        if not lines or not lines[0].strip():
            _warnEmptyOutput(cmd)
            content = self._FALLBACK
        return _GuardedPipe(content)


def _warnEmptyOutput(cmd):
    """ Capped per process (each forked worker gets its own counter): if something
    systematic is broken - e.g. a Vina build that refuses to --score_only without a grid
    box - this fires once per conformer, and an uncapped log would be tens of MB. """
    global _emptyOutputWarnings
    _emptyOutputWarnings += 1
    if _emptyOutputWarnings > _MAX_EMPTY_OUTPUT_WARNINGS:
        return
    message = '[protac] WARNING: no output from shell command, using 0 instead: %s\n' % cmd
    if _emptyOutputWarnings == _MAX_EMPTY_OUTPUT_WARNINGS:
        message += ('[protac] WARNING: further empty-output warnings from this process '
                    'are suppressed.\n')
    _warn(message)


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
