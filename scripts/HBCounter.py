# -*- coding: utf-8 -*-


# Standard imports
import os
import argparse as ap
from multiprocessing import Pool
from functools import partial
from pathlib import Path
from typing import List, Tuple, Set, Dict

# Local imports
from Helpers import SimIt, HBondLinker
from Helpers import hbond_mod as hbm
from Helpers.ReportUtils import extract_metrics

# External imports
import mdtraj as md
import pandas as pd
import numpy as np


# Script information
__author__ = "Marti Municoy, Carles Perez"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy, Carles Perez"
__email__ = "marti.municoy@bsc.es, carles.perez@bsc.es"


def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument("sim_paths", metavar="PATH", type=str,
                        nargs='*',
                        help="Path to PELE simulation folders")

    parser.add_argument("-l", "--ligand_resname",
                        metavar="LIG", type=str, default='LIG',
                        help="Ligand residue name")

    parser.add_argument("-d", "--distance",
                        metavar="D", type=float, default='0.25',
                        help="Hydrogen bonds distance")

    parser.add_argument("-a", "--angle",
                        metavar="A", type=float, default='2.0943951023931953',
                        help="Hydrogen bonds angle")

    parser.add_argument("-p", "--pseudo_hb",
                        metavar="BOOL", type=bool, default=False,
                        help="Look for pseudo hydrogen bonds")

    parser.add_argument("-n", "--processors_number",
                        metavar="N", type=int, default=None,
                        help="Number of processors")

    parser.add_argument("-t", "--topology_path",
                        metavar="PATH", type=str,
                        default='output/topologies/topology_0.pdb',
                        help="Relative path to topology")

    parser.add_argument("-o", "--output_path",
                        metavar="PATH", type=str,
                        default='hbonds.csv',
                        help="Relative path to output file")

    parser.add_argument("-pdb", "--pdb_file",
                        metavar="BOOL", type=bool,
                        default=False,
                        help="Set if the output is in PDB format.")

    parser.add_argument("--PELE_output_path",
                        metavar="PATH", type=str, default='output',
                        help="Relative path to PELE output folder")

    parser.add_argument("-r", "--report_name",
                        metavar="NAME", type=str,
                        default='report',
                        help="PELE report name")

    parser.add_argument('--include_rejected_steps',
                        dest='include_rejected_steps',
                        action='store_true')

    parser.add_argument("--alternative_output_path",
                        metavar="PATH", type=str, default=None,
                        help="Alternative path to save output results")
   
    parser.add_argument("-ad", "--adaptive",
                        metavar="BOOL", type=bool,
                        default=True,
                        help="Set to False if adaptive is not set.")

    parser.set_defaults(include_rejected_steps=False)

    args = parser.parse_args()

    return args.sim_paths, args.ligand_resname, args.distance, args.angle, \
        args.pseudo_hb, args.processors_number, args.topology_path, \
        args.output_path, args.report_name, args.PELE_output_path, \
        args.include_rejected_steps, args.alternative_output_path, \
        args.pdb_file, args.adaptive


def account_for_ignored_hbonds(hbonds_in_traj: List[List[HBondLinker]],
                               accepted_steps: List[int]
                               ) -> List[List[HBondLinker]]:
    new_hbonds_in_traj = []
    for i, s in enumerate(accepted_steps):
        new_hbonds_in_traj.append(hbonds_in_traj[s])

    return new_hbonds_in_traj


def find_hbonds_in_trajectory(lig_resname: str, distance: float, angle: float,
                              pseudo: bool, topology_path: Path,
                              chain_ids: List[str], report_name: str,
                              include_rejected_steps: bool, traj_path: Path
                              ) -> Tuple[List[List[HBondLinker]],
                                         List[int], List[int],
                                         Set[md.core.topology.Atom],
                                         Set[md.core.topology.Atom]]:
    if os.path.splitext(str(traj_path))[-1] == '.pdb':
        pdb = True
    else:
        pdb = False
    try:
        if not pdb:
            traj = md.load_xtc(str(traj_path), top=str(topology_path))
        else:
            traj = md.load_pdb(str(traj_path))
    except OSError:
        print('     - Warning: problems loading trajectory '
              '{}, it will be ignored'.format(traj_path))
        # Return empty data structures
        return list(), list(), list(), set(), set()
    lig = traj.topology.select('resname {}'.format(lig_resname))
    hbonds_in_traj, donors, acceptors = find_ligand_hbonds(traj, lig, distance,
                                                           angle, pseudo,
                                                           chain_ids)

    # Recover corresponding report
    num = int(''.join(filter(str.isdigit, traj_path.name)))
    print(num)
    path = traj_path.parent
    print(path)
    report_path = path.joinpath(report_name + '_{}'.format(num))
    print(str(report_path))
    metrics = extract_metrics((report_path, ), (2, 3))[0]

    total_steps = []
    accepted_steps = []
    for m in metrics:
        total_steps.append(int(m[0]))
        accepted_steps.append(int(m[1]))

    try:
        if (include_rejected_steps):
            hbonds_in_traj = account_for_ignored_hbonds(hbonds_in_traj,
                                                        accepted_steps)
        else:
            if (len(hbonds_in_traj) != len(total_steps)
                    or len(hbonds_in_traj) != len(accepted_steps)):
                raise IndexError
    except IndexError:
        print('     - Warning: inconsistent number of models found in '
              + 'trajectory {} from {}, '.format(num, path)
              + 'this trajectory will be ignored')
        # Return empty data structures
        return list(), list(), list(), set(), set()

    return hbonds_in_traj, total_steps, accepted_steps, donors, acceptors


def find_ligand_hbonds(traj: md.Trajectory, lig: np.ndarray, distance: float,
                       angle: float, pseudo: bool, chain_ids: List[str]
                       ) -> Tuple[List[List[HBondLinker]],
                                  Set[md.core.topology.Atom],
                                  Set[md.core.topology.Atom]]:
    hbonds_list = []
    donors = set()
    acceptors = set()
    for model_id in range(0, traj.n_frames):
        results, _donors, _acceptors = find_hbond_in_snapshot(
            traj, model_id, lig, distance, angle, pseudo, chain_ids)
        hbonds_list.append(results)
        for d in _donors:
            donors.add(d)
        for a in _acceptors:
            acceptors.add(a)

    return hbonds_list, donors, acceptors


def find_hbond_in_snapshot(traj: md.Trajectory, model_id: int, lig: np.ndarray,
                           distance: float, angle: float, pseudo: bool,
                           chain_ids: List[str]
                           ) -> Tuple[List[HBondLinker],
                                      Set[md.core.topology.Atom],
                                      Set[md.core.topology.Atom]]:
    hbonds = hbm.baker_hubbard(traj=traj[model_id], distance=distance,
                               angle=angle, pseudo=pseudo)

    results = []
    donors = set()
    acceptors = set()
    for hbond in hbonds:
        if (hbond[0] in lig):
            donors.add(traj.topology.atom(hbond[0]))
        elif (hbond[2] in lig):
            acceptors.add(traj.topology.atom(hbond[2]))
        if (any(atom in lig for atom in hbond) and not
                all(atom in hbond for atom in lig)):
            for atom in hbond:
                if (atom not in lig):
                    _atom = traj.topology.atom(atom)
                    hb_linker = HBondLinker(
                        chain_ids[_atom.residue.chain.index],
                        _atom.residue, tuple((_atom.name, )))
                    results.append(hb_linker)
                    break

    return results, donors, acceptors


def parse_results(results: Tuple[List[HBondLinker],
                                 Set[md.core.topology.Atom],
                                 Set[md.core.topology.Atom]],
                  trajectories: List[md.Trajectory]
                  ) -> Tuple[pd.DataFrame,
                             Set[md.core.topology.Atom],
                             Set[md.core.topology.Atom],
                             int]:
    counter = 0
    data = pd.DataFrame()
    donors = set()
    acceptors = set()
    for (r, t_steps, a_steps, _donors, _acceptors), t in zip(results,
                                                             trajectories):
        counter += len(r)
        try:
            epoch = int(t.parent.name)
        except ValueError:
            epoch = 0 
        trajectory = int(''.join(filter(str.isdigit, t.name)))
        for hbonds, t_s, a_s in zip(r, t_steps, a_steps):
            data = data.append(pd.DataFrame(
                [(epoch, trajectory, t_s, a_s, hbonds)],
                columns=['epoch', 'trajectory', 'step', 'model',
                         'hbonds']))

        for d in _donors:
            donors.add(d)
        for a in _acceptors:
            acceptors.add(a)

    return data, donors, acceptors, counter


def main():
    # Parse args
    PELE_sim_paths, lig_resname, distance, angle, pseudo_hb, proc_number, \
        topology_relative_path, output_relative_path, report_name, \
        PELE_output_path, include_rejected_steps, alternative_output_path, \
        pdb_file, adaptive = parse_args()

    all_sim_it = SimIt(PELE_sim_paths)
    print(' - The following PELE simulation paths will be analyzed:')
    print(' Adaptive flag: {}'.format(adaptive))
    for PELE_sim_path in all_sim_it:
        print('   - {}'.format(PELE_sim_path))

    for PELE_sim_path in all_sim_it:
        print(' - Analyzing {}'.format(PELE_sim_path))
        topology_path = PELE_sim_path.joinpath(topology_relative_path)
        print(topology_path)
        if (not topology_path.is_file()):
            print(' - Skipping simulation because topology file with '
                  + 'connectivity was missing')
            continue

        # Retrieve chain ids
        chain_ids = set()
        with open(str(topology_path), 'r') as f:
            for line in f:
                if (len(line) < 80):
                    continue
                line = line.strip()
                chain_ids.add(line[21])
        chain_ids = sorted(list(chain_ids))

        parallel_function = partial(find_hbonds_in_trajectory, lig_resname,
                                    distance, angle, pseudo_hb, topology_path,
                                    chain_ids, report_name,
                                    include_rejected_steps)

        sim_it = SimIt(PELE_sim_path)
        if not pdb_file:
            sim_it.build_traj_it(PELE_output_path, 'trajectory', 'xtc', adaptive=adaptive)
        else:
            sim_it.build_traj_it(PELE_output_path, 'trajectory', 'pdb', adaptive=adaptive)

        trajectories = [traj for traj in sim_it.traj_it]
        print(trajectories)
        with Pool(proc_number) as pool:
            results = pool.map(parallel_function,
                               trajectories)

        data, donors, acceptors, counter = parse_results(results, trajectories)

        print('     - {} models were found'.format(counter))
        print('     - {} ligand donors were found'.format(len(donors)))
        print('     - {} ligand acceptors were found'.format(len(acceptors)))

        if (alternative_output_path is not None):
            output_path = Path(alternative_output_path)
            output_path = output_path.joinpath(PELE_sim_path.name)
            output_path = output_path.joinpath(output_relative_path)
            try:
                os.makedirs(str(output_path.parent))
            except FileExistsError:
                pass
        else:
            output_path = PELE_sim_path.joinpath(output_relative_path)

        output_info_path = Path(output_path.parent).joinpath(
            output_path.name.replace(output_path.suffix, '') + '.info')

        with open(str(output_info_path), 'w') as file:
            file.write(str(PELE_sim_path.name) + '\n')
            file.write('{} donors: {}\n'.format(len(donors),
                                                list(donors)))
            file.write('{} acceptors: {}\n'.format(len(acceptors),
                                                   list(acceptors)))

        data.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
