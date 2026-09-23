#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Python 3 driver for PRosettaC (LondonLab/PRosettaC). Launched by
# protac/protocols/protocol_prosettac.py under PROSETTAC_PYTHON_DIC's conda env, always
# from the work dir (extra/prosettac/) - every phase reads/writes fixed relative
# filenames there, replicating main.py's own layout. Only stages inputs and calls
# PRosettaC's own utils/rosetta/protac_lib functions - none of their logic is
# reimplemented.
#
# Covers all 7 phases of the original main.py. The two phases the original farms out to a
# PBS/SGE/SLURM scheduler (local docking, constrained linker conformations) run here in a
# local thread pool instead - see _runCommandsInParallel().

import argparse
import concurrent.futures
import glob
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.environ['SCRIPTS_FOL'].rstrip(os.sep))

STATE_FILE = 'state.json'
PATCHDOCK_RESULTS = 'Patchdock_Results'
PROTAC_SMI = 'protac.smi'
# Lines of a failed job's output to quote when reporting it.
FAILURE_TAIL_LINES = 30


def _readState():
    with open(STATE_FILE) as f:
        return json.load(f)


def _writeState(state):
    with open(STATE_FILE, 'w') as f:
        json.dump(state, f)


def _stageSimpleName(srcPath, destName):
    """ Copies srcPath to destName in cwd. Several PRosettaC functions derive filenames
    with "<path>.split('.')[0]", which breaks on any dot outside the extension - staging
    under a dot-free name avoids that regardless of where the original file lives. """
    shutil.copy(srcPath, destName)
    return destName


def _normalizeHead(sdfFile, anchor):
    """ Rewrites sdfFile as heavy atoms only, no radicals, aromatic bonds; returns the
    0-based anchor re-indexed accordingly. translate_anchors() compares the file read
    unsanitized against its sanitized H-added copy, so explicit H, Kekule rings or the
    radicals of a .pdbqt-derived SDF all make it return -1. """
    from rdkit import Chem

    mol = Chem.SDMolSupplier(sdfFile, removeHs=False, sanitize=False)[0]
    if mol is None:
        raise RuntimeError(f'RDKit cannot read the warhead file {sdfFile}.')
    if not 0 <= anchor < mol.GetNumAtoms():
        raise RuntimeError(f'Anchor atom {anchor + 1} is out of range: the warhead has '
                           f'{mol.GetNumAtoms()} atoms.')
    if mol.GetAtomWithIdx(anchor).GetAtomicNum() == 1:
        raise RuntimeError(f'Anchor atom {anchor + 1} is a hydrogen; pick a heavy atom.')
    newAnchor = sum(1 for a in mol.GetAtoms() if a.GetIdx() < anchor and a.GetAtomicNum() != 1)
    for atom in mol.GetAtoms():
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(False)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    Chem.MolToMolFile(mol, sdfFile, kekulize=False)
    return newAnchor


def _catFiles(parts, dest):
    with open(dest, 'wb') as out:
        for part in parts:
            with open(part, 'rb') as f:
                shutil.copyfileobj(f, out)
    return dest


def _removeStale(directory, *patterns):
    """ Deletes what a previous run of the same phase left behind. Not cosmetic: the
    original refuses to overwrite in places (mol_to_params without --clobber) and appends
    in others (Rosetta score files), so a relaunched Scipion step would otherwise either
    die or mix old and new models. """
    for pattern in patterns:
        for path in glob.glob(os.path.join(directory, pattern)):
            os.remove(path)


def _runShellCommand(command, cwd):
    proc = subprocess.run(command, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, universal_newlines=True)
    return proc.returncode, proc.stdout


def _runCommandsInParallel(commands, cwd, threads, what):
    """ Local stand-in for the pair of calls main.py makes for phases 5 and 6,
    'cluster.runBatchCommands(commands, mem)' + 'cluster.wait(jobs)'. The commands are
    the ones the original builds, unchanged - only what runs them differs (a bounded
    thread pool here, a PBS/SGE/SLURM submission there). No PRosettaC logic is bypassed
    by this: rs.local_docking() merely returns a command string and
    constraint_generation.py is a standalone script; neither ever touches the cluster
    layer, only main.py does.

    A failing job is reported but doesn't abort the phase, matching the original, where
    each command is an independent cluster job and the later phases just work with
    whatever models exist. All of them failing does raise: nothing downstream can work. """
    failures = []
    total = len(commands)
    print('[protac] %s: %d job(s), %d at a time' % (what, total, threads))
    sys.stdout.flush()

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        futures = {pool.submit(_runShellCommand, cmd, cwd): cmd for cmd in commands}
        for done, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            returnCode, output = future.result()
            if returnCode != 0:
                failures.append((futures[future], output))
            if done % 50 == 0 or done == total:
                print('[protac] %s: %d/%d done, %d failed'
                      % (what, done, total, len(failures)))
                sys.stdout.flush()

    if failures:
        command, output = failures[0]
        message = ('%d of %d %s job(s) failed. First failure:\n  %s\n%s'
                   % (len(failures), total, what, command,
                      '\n'.join(output.splitlines()[-FAILURE_TAIL_LINES:])))
        if len(failures) == total:
            raise RuntimeError(message)
        print('[protac] WARNING: ' + message)
        sys.stdout.flush()


def runPrepare(args):
    """ Phases 1+2 (entangled in main.py, not separable): per-structure addH_sdf ->
    translate_anchors -> clean -> mol_to_params -> relax -> clean. """
    import rosetta as rs
    import protac_lib as pl
    import utils

    structs = [_stageSimpleName(args.struct1, 'Struct0.pdb'),
              _stageSimpleName(args.struct2, 'Struct1.pdb')]
    chains = [args.chain1, args.chain2]
    heads = [_stageSimpleName(args.head1, 'Head0.sdf'),
            _stageSimpleName(args.head2, 'Head1.sdf')]
    anchors = [args.anchor1 - 1, args.anchor2 - 1]  # form is 1-based, PRosettaC is 0-based

    ptParams = []
    for i in (0, 1):
        newHead = f'Head{i}_H.sdf'
        anchors[i] = _normalizeHead(heads[i], anchors[i])
        utils.addH_sdf(heads[i], newHead)
        anchors[i] = pl.translate_anchors(heads[i], newHead, anchors[i])
        if anchors[i] == -1:
            raise RuntimeError(
                f'translate_anchors failed for head {i + 1}: OpenBabel re-perceived the '
                'warhead as a different molecule. Check that the warhead SDF is 3D '
                '(a 2D file loses its bond orders on the way through PDB).')
        heads[i] = newHead

        rs.clean(structs[i], chains[i])
        structs[i] = f"{structs[i].split('.')[0]}_{chains[i]}.pdb"
        ptPdb, ptParam = rs.mol_to_params(heads[i], f'PT{i}', f'PT{i}')
        ptParams.append(ptParam)
        _catFiles([ptPdb, structs[i]], f'Side{i}.pdb')
        rs.relax(f'Side{i}.pdb', ptParam)
        rs.clean(f'Side{i}_0001.pdb', chains[i])
        structs[i] = f'Init{i}.pdb'
        _catFiles([ptPdb, f'Side{i}_0001_{chains[i]}.pdb'], structs[i])

    _writeState({'heads': heads, 'anchors': anchors, 'initFiles': structs,
                'ptParams': ptParams})


def runSampleDist(args):
    """ Phase 3: linker distance sampling, reads runPrepare's own state.json. """
    import protac_lib as pl

    state = _readState()
    with open(PROTAC_SMI, 'w') as f:
        f.write(args.smiles.strip() + '\n')

    minValue, maxValue = pl.SampleDist(state['heads'], state['anchors'], PROTAC_SMI)
    if (minValue, maxValue) == (None, None):
        raise RuntimeError(
            'SampleDist found no substructure match between a head .sdf and the PROTAC '
            'SMILES. Check that both warhead SDFs are exact substructures of "PROTAC '
            'SMILES", in their bound conformation.')
    if (minValue, maxValue) == (0, 0):
        raise RuntimeError(
            'SampleDist could not generate any PROTAC conformation to sample the anchor '
            'distance. Check that both warhead SDFs are in a valid bound conformation and '
            'that "Structure 1/2 anchor atom" point at the right atoms.')

    state.update({'minValue': minValue, 'maxValue': maxValue})
    _writeState(state)


def runPatchdock(args):
    """ Phase 4: PatchDock global docking, reads runSampleDist's own state.json. """
    import utils

    state = _readState()
    oneBasedAnchors = [a + 1 for a in state['anchors']]
    numResults = utils.patchdock(state['initFiles'], oneBasedAnchors,
                                 state['minValue'], state['maxValue'],
                                 args.globalResults, args.threshold)
    if numResults is None:
        raise RuntimeError(
            'PatchDock found no global docking solution within the geometrical '
            'constraints (the anchor distance window from the sampling step).')

    state['numResults'] = numResults
    _writeState(state)


def runLocalDocking(args):
    """ Phase 5: RosettaScripts local docking of every PatchDock result, one job per
    result. Each command line comes from PRosettaC's own rs.local_docking(), which is all
    that function does - main.py is what submits them, and _runCommandsInParallel() is
    what takes that role here. """
    import rosetta as rs

    state = _readState()
    # The original passes these as absolute paths too (main.py: curr_dir + '/' + ...):
    # the jobs run from Patchdock_Results/, the .params files live one level up.
    ptParams = [os.path.abspath(p) for p in state['ptParams']]

    # local.fasc is appended to by every job, and a rerun that yields fewer solutions
    # would leave stale *_docking_????.pdb behind for the next phase to pick up.
    _removeStale(PATCHDOCK_RESULTS, 'local.fasc', '*_docking_????.pdb')

    commands = [rs.local_docking('pd.%d.pdb' % (i + 1), args.chain1 + 'X',
                                 args.chain2 + 'Y', ptParams[0], ptParams[1],
                                 args.nstruct)
               for i in range(state['numResults'])]
    _runCommandsInParallel(commands, PATCHDOCK_RESULTS, args.threads, 'local docking')


def _dockingSuffix(solution):
    """ 'pd.<i>_docking_<nnnn>.pdb' -> '<i>_<n>', the per-job output suffix main.py
    derives with "s.split('.')[1].split('_')". Reproduced exactly: constraint_generation.py
    names every file it writes after it, and clustering.py rebuilds the same two numbers
    back out of those names. """
    parts = solution.split('.')[1].split('_')
    return '%s_%d' % (parts[0], int(parts[2]))


def runConstraintConf(args):
    """ Phase 6: up to 100 constrained conformations of the full linker per local docking
    solution, by running PRosettaC's own constraint_generation.py once per solution - the
    same commands main.py builds, run the same way as phase 5. """
    state = _readState()
    scriptsFol = os.environ['SCRIPTS_FOL'].rstrip(os.sep)

    solutions = sorted(os.path.basename(p) for p in
                      glob.glob(os.path.join(PATCHDOCK_RESULTS, '*_docking_????.pdb')))
    if not solutions:
        raise RuntimeError(
            'No local docking solution (%s/*_docking_????.pdb) to build linker '
            'conformations from: the local docking phase produced nothing usable. Check '
            'its log for Rosetta errors.' % PATCHDOCK_RESULTS)

    # Everything constraint_generation.py writes, cleared for a rerun: it calls
    # mol_to_params with overwrite=False (which then refuses to overwrite its output) and
    # the relax at its end appends to score.sc.
    _removeStale(PATCHDOCK_RESULTS, 'score.sc', 'confs_*.sdf', 'v_*.sdf', 'docked_*.pdb',
                 'docked_*.sdf', 'PT_*.pdb', 'PT_*.params', 'combined_*.pdb')

    chains = args.chain1 + args.chain2
    script = os.path.join(scriptsFol, 'constraint_generation.py')
    # sys.executable, not the bare 'python' main.py uses: PATH carries both this env's
    # bin (prepended by the conda activation that launched this driver) and
    # PROSETTAC_PYTHON2_HOME/bin (prepended by getPRosettaCEnviron, for mol_to_params'
    # bare 'python2.7'), so which 'python' wins depends on the order those two are
    # applied in. constraint_generation.py needs this Python 3 + RDKit one; naming the
    # interpreter outright takes the question off the table.
    commands = ['%s %s ../%s ../%s ../%s %s %s %s'
                % (sys.executable, script, state['heads'][0], state['heads'][1],
                   PROTAC_SMI, _dockingSuffix(s), s, chains)
               for s in solutions]
    _runCommandsInParallel(commands, PATCHDOCK_RESULTS, args.threads,
                           'constrained conformation generation')


def runClustering(args):
    """ Phase 7: DBSCAN clustering of the top models by interface RMSD, calling
    PRosettaC's own clustering.main() with the 4 arguments main.py passes it on the
    command line. """
    import clustering

    # clustering.py compares the PDB chain column one character at a time (line[21] ==
    # chain), so a multi-chain moving body matches no CA atom at all and it dies on a
    # division by zero instead of saying what is wrong. _validate() rejects this too.
    if len(args.chain2) != 1:
        raise RuntimeError(
            'The moving structure must have a single chain ID for clustering, got "%s". '
            "PRosettaC's clustering.py matches the PDB chain column against it character "
            'by character.' % args.chain2)

    _catFiles(['Init0.pdb', 'Init1.pdb'], 'Init.pdb')

    workDir = os.getcwd()
    try:
        # clustering.main() chdir's into Patchdock_Results and back, but returns early -
        # still inside it - when there is nothing to cluster, hence the restore.
        clustering.main('clustering.py', [str(args.topScore), str(args.topLocal),
                                          str(args.rmsd), args.chain2])
    finally:
        os.chdir(workDir)

    if not os.path.isdir('Results'):
        raise RuntimeError(
            'Clustering produced no Results/ directory. Either %s/score.sc was never '
            'written (no local docking solution survived the linker conformation phase, '
            'i.e. no valid PROTAC conformation could bridge any of them), or no model '
            'scored below the 0 energy threshold.' % PATCHDOCK_RESULTS)


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', required=True,
                        choices=['prepare', 'sampledist', 'patchdock', 'localdocking',
                                 'constraintconf', 'clustering'])
    parser.add_argument('--struct1'); parser.add_argument('--chain1')
    parser.add_argument('--struct2'); parser.add_argument('--chain2')
    parser.add_argument('--head1'); parser.add_argument('--anchor1', type=int)
    parser.add_argument('--head2'); parser.add_argument('--anchor2', type=int)
    parser.add_argument('--smiles')
    parser.add_argument('--global-results', dest='globalResults', type=int, default=500)
    parser.add_argument('--threshold', type=float, default=2.0)
    parser.add_argument('--nstruct', type=int, default=10)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--top-score', dest='topScore', type=int, default=1000)
    parser.add_argument('--top-local', dest='topLocal', type=int, default=200)
    parser.add_argument('--cluster-rmsd', dest='rmsd', type=float, default=4.0)
    return parser.parse_args()


if __name__ == '__main__':
    args = parseArgs()
    if args.phase == 'prepare':
        runPrepare(args)
    elif args.phase == 'sampledist':
        runSampleDist(args)
    elif args.phase == 'patchdock':
        runPatchdock(args)
    elif args.phase == 'localdocking':
        runLocalDocking(args)
    elif args.phase == 'constraintconf':
        runConstraintConf(args)
    elif args.phase == 'clustering':
        runClustering(args)
