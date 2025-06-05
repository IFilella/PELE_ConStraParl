# -*- coding: utf-8 -*-


# Standard imports
import os
import argparse as ap
from multiprocessing import Pool
from functools import partial
from collections import defaultdict
from operator import itemgetter

# Local imports
from Helpers import SimIt
from Helpers.ReportUtils import extract_PELE_ids
from Helpers.ReportUtils import get_metric_by_PELE_id
from Helpers.Hbonds import (
    extract_hbond_linkers, hbond_fulfillment, get_hbond_linkers,
    check_hbonds_linkers, print_hbonds
)

# External imports
import mdtraj as md
import numpy as np
from sklearn.cluster import MeanShift
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt


# Script information
__author__ = "Marti Municoy, Carles Perez"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy, Carles Perez"
__email__ = "marti.municoy@bsc.es, carles.perez@bsc.es"


def parse_args():
    parser = ap.ArgumentParser()
    parser.add_argument("traj_paths", metavar="PATH", type=str,
                        nargs='*',
                        help="Path to PELE trajectory files")

    parser.add_argument("-l", "--ligand_resname",
                        metavar="LIG", type=str, default='LIG',
                        help="Ligand residue name")

    parser.add_argument("-b", "--bandwidth",
                        metavar="B", type=float, default='5',
                        help="Clustering bandwidth")

    parser.add_argument("-n", "--processors_number",
                        metavar="N", type=int, default=None,
                        help="Number of processors")

    parser.add_argument("--ie_col",
                        metavar="N", type=int, default=5,
                        help="Interaction energy column in PELE reports")

    parser.add_argument("--rmsd_col",
                        metavar="N", type=int, default=7,
                        help="RMSD column in PELE reports")

    parser.add_argument("-t", "--topology_path",
                        metavar="PATH", type=str,
                        default='output/topologies/topology_0.pdb',
                        help="Relative path to topology")

    parser.add_argument("-r", "--report_name",
                        metavar="NAME", type=str,
                        default='report',
                        help="PELE report name")

    parser.add_argument("--hbonds_path",
                        metavar="PATH", type=str,
                        default='hbonds.csv',
                        help="Path to H bonds file")

    parser.add_argument("-g1", "--golden_hbonds_1", nargs='*',
                        metavar="C:R[:A1, A2]", type=str, default=[],
                        help="Chain (C), residue (R) [and atoms (A1, A2)] of"
                        + "subset 1 of golden H bonds. Subset 1 contains H "
                        + "bond conditions that must always be fulfilled in "
                        + "the filtering process")

    parser.add_argument("-g2", "--golden_hbonds_2", nargs='*',
                        metavar="C:R[:A1, A2]", type=str, default=[],
                        help="Chain (C), residue (R) [and atoms (A1, A2)] of"
                        + "subset 2 of golden H bonds. Subset 2 contains H "
                        + "bond conditions that only a minimum number of them "
                        + "must be fulfilled in the filtering process. The "
                        + "minimum of required conditions from subset 2 is "
                        "defined with the minimum_g2_conditions argument")

    parser.add_argument("--minimum_g2_conditions",
                        metavar="N", type=int, default=2,
                        help="Minimum number of subset 2 golden H bonds "
                        + "that must be fulfilled in the filtering process")

    parser.add_argument("-o", "--output_path",
                        metavar="PATH", type=str, default="filtering_results")

    parser.add_argument('--generate_plots', dest='generate_plots',
                        action='store_true')

    parser.add_argument('--representative_extraction_method',
                        choices=["interaction_energy",
                                 "center",
                                 "total_energy"],
                        type=str, metavar="MODE",
                        default="interaction_energy",
                        help="Method to extract the representative structure"
                        + " from the final selected cluster. The structure can"
                        + " be selected by minimizing the interaction energy, "
                        + "the total energy or by considering the one that is "
                        + "closer to the center of the cluster.")

    parser.add_argument("--maximum_cluster_density_threshold",
                        metavar="B", type=float, default='0.25',
                        help="Maximum cluster density threshold that is used"
                        + " to filter out low populated clusters in the final"
                        + " selection stage.")

    parser.set_defaults(generate_plots=False)

    args = parser.parse_args()

    return args.traj_paths, args.ligand_resname, args.bandwidth, \
        args.processors_number, args.ie_col, args.rmsd_col, \
        args.topology_path, args.report_name, args.hbonds_path, \
        args.golden_hbonds_1, args.golden_hbonds_2, \
        args.minimum_g2_conditions, args.output_path, args.generate_plots, \
        args.representative_extraction_method, \
        args.maximum_cluster_density_threshold


def get_reports_list(trajectories, report_name):
    reports_list = []

    for traj in trajectories:
        # Recover corresponding report
        num = int(''.join(filter(str.isdigit, traj.name)))
        path = traj.parent

        report_path = path.joinpath(report_name + '_{}'.format(num))

        reports_list.append(report_path)

    return reports_list


def extract_ligand_properties(topology_path, resname):
    n_heavy_atoms = 0
    molecular_weight = 0

    traj = md.load_pdb(str(topology_path))
    lig = traj.topology.select('resname {}'.format(resname))
    for atom in traj.atom_slice(lig).top.atoms:
        if (str(atom.element) != 'hydrogen'):
            n_heavy_atoms += 1
        molecular_weight += atom.element.mass

    return n_heavy_atoms, molecular_weight


def p_extract_ligand_metrics(cols, report_path):
    results = []

    if (report_path.is_file()):
        try:
            with open(str(report_path), 'r') as f:
                f.readline()
                for line in f:
                    line.strip()
                    fields = line.split()
                    metrics = []
                    for col in cols:
                        metrics.append(fields[col - 1])
                    results.append(metrics)

        except IndexError:
            print(' - p_extract_ligand_metric Warning: wrong index '
                  + 'supplied for trajectory: \'{}\''.format(report_path))
    else:
        print(' - p_extract_ligand_metric Warning: wrong path to report '
              + 'for trajectory: \'{}\''.format(report_path))

    return results


def extract_ligand_metrics(reports, cols, proc_number):
    parallel_function = partial(p_extract_ligand_metrics, cols)

    with Pool(proc_number) as pool:
        results = pool.map(parallel_function,
                           reports)

    return results


def filter_by_hbonds(f_hbond_data):
    # Compatibility with older pandas versions
    try:
        unhashable_PELE_ids = f_hbond_data.loc[:, ('epoch', 'trajectory',
                                                   'model')].to_list()
    except AttributeError:
        unhashable_PELE_ids = f_hbond_data.loc[:, ('epoch', 'trajectory',
                                                   'model')].values

    PELE_ids = []
    for PELE_id in unhashable_PELE_ids:
        PELE_ids.append(tuple(PELE_id))

    return PELE_ids


def filter_by_energies(PELE_ids, ie_by_PELE_id):
    ies = []
    for PELE_id in PELE_ids:
        ies.append(ie_by_PELE_id[PELE_id])

    upper_bound = np.percentile(ies, 25)
    print('   - Energetic upper bound: {:.1f} kcal/mol'.format(upper_bound))

    filtered_PELE_ids = []
    for PELE_id in PELE_ids:
        ie = ie_by_PELE_id[PELE_id]

        if (ie < upper_bound):
            filtered_PELE_ids.append(PELE_id)

    return filtered_PELE_ids


def extract_ligand_coords(filtered_PELE_ids, trajectories, lig_resname,
                          topology_path, proc_number):
    PELE_ids_dict = defaultdict(list)
    for (e, t, m) in filtered_PELE_ids:
        PELE_ids_dict[(e, t)].append(m)

    filtered_traj = []
    for traj in trajectories:
        e = int(traj.parent.name)
        t = int(''.join(filter(str.isdigit, traj.name)))
        if ((e, t) in PELE_ids_dict):
            filtered_traj.append(traj)

    parallel_function = partial(p_extract_ligand_coords,
                                lig_resname,
                                topology_path)

    with Pool(proc_number) as pool:
        results = pool.map(parallel_function,
                           filtered_traj)

    filtered_results = []
    filtered_PELE_ids = []
    for i, traj in enumerate(filtered_traj):
        e = int(traj.parent.name)
        t = int(''.join(filter(str.isdigit, traj.name)))
        for m in PELE_ids_dict[(e, t)]:
            filtered_results.append(results[i][m])
            filtered_PELE_ids.append((e, t, m))

    return filtered_results, filtered_PELE_ids


def p_extract_ligand_coords(lig_resname, topology_path, traj_path):
    try:
        traj = md.load_xtc(str(traj_path), top=str(topology_path))
    except OSError:
        print('     - Warning: problems loading trajectory '
              '{}, it will be ignored'.format(traj_path))
        return {}

    lig = traj.topology.select('resname {}'.format(lig_resname))

    lig_coords = traj.atom_slice(lig).xyz

    reshaped_lig_coords = []
    for chunk in lig_coords:
        reshaped_lig_coords.append(np.concatenate(chunk))

    return reshaped_lig_coords


def clusterize(lig_coords, bandwidth, proc_number):
    gm = MeanShift(bandwidth=bandwidth,
                   n_jobs=proc_number,
                   cluster_all=True)

    return gm.fit_predict(lig_coords), gm.cluster_centers_


def calculate_probabilities(cluster_results):
    p_dict = defaultdict(int)
    total = 0
    for cluster in cluster_results:
        p_dict[cluster] += 1
        total += 1

    for cluster, probability in p_dict.items():
        p_dict[cluster] /= total

    return p_dict


def get_most_populated_clusters(p_dict, cluster_centers_array,
                                maximum_cluster_density_threshold):
    ordered_clusters = sorted(p_dict.items(), key=itemgetter(1),
                              reverse=True)

    # Prevent filtering all clusters out in case were small densities are found
    population_threshold = min((maximum_cluster_density_threshold,
                                ordered_clusters[0][1]))

    print(' - Selecting clusters with density higher than {:.2f}'.format(
        population_threshold))

    # Allow to count all clusters if there is a draw
    selected_clusters = []
    for c, p in sorted(p_dict.items(), key=itemgetter(1),
                       reverse=True):
        if (p >= population_threshold):
            selected_clusters.append(c)

    selected_centers = []
    for c in selected_clusters:
        selected_centers.append(cluster_centers_array[c])

    print('   - Number of selected clusters: {}'.format(
        len(selected_clusters)))

    return selected_clusters, selected_centers


def get_best_cluster(cluster_ids, selected_centers, results, filtered_PELE_ids,
                     ie_by_PELE_id, rmsd_by_PELE_id, te_by_PELE_id, lig_coords,
                     representative_extraction_method):
    cluster_metrics = {}
    for cluster_id, center in zip(cluster_ids, selected_centers):
        ies, rmsds, tes, representative_PELE_id = \
            get_metrics_from_cluster(cluster_id, results,
                                     filtered_PELE_ids, ie_by_PELE_id,
                                     rmsd_by_PELE_id, te_by_PELE_id,
                                     lig_coords, center,
                                     representative_extraction_method)
        cluster_metrics[cluster_id] = (ies, rmsds, tes, representative_PELE_id)

    if (len(cluster_ids) == 1):
        return cluster_id, ies, rmsds, tes, representative_PELE_id

    mean_ie = np.mean(ies)
    for current_cluster_id, metrics in cluster_metrics.items():
        if (np.mean(metrics[0]) < mean_ie):
            cluster_id = current_cluster_id
            ies, rmsds, tes, representative_PELE_id = metrics

    print('   - Selected cluster {} (lowest mean interaction energy)'.format(
        cluster_id))

    return cluster_id, ies, rmsds, tes, representative_PELE_id


def get_metrics_from_cluster(cluster_id, results, PELE_ids,
                             ie_by_PELE_id, rmsd_by_PELE_id,
                             te_by_PELE_id, lig_coords,
                             cluster_center,
                             representative_extraction_method):
    ies = []
    rmsds = []
    tes = []
    best_metric = None
    best_id = None
    for i, (r, PELE_id) in enumerate(zip(results, PELE_ids)):
        if (r == cluster_id):
            ies.append(ie_by_PELE_id[PELE_id])
            rmsds.append(rmsd_by_PELE_id[PELE_id])
            tes.append(te_by_PELE_id[PELE_id])

            if (representative_extraction_method == "center"):
                diff = lig_coords[i] - cluster_center
                sq_dist = np.dot(diff, diff)
                if (best_metric is None or sq_dist < best_metric):
                    best_metric = sq_dist
                    best_id = PELE_id

            elif (representative_extraction_method == "interaction_energy"):
                if (best_metric is None
                        or ie_by_PELE_id[PELE_id] < best_metric):
                    best_metric = ie_by_PELE_id[PELE_id]
                    best_id = PELE_id

            elif (representative_extraction_method == "total_energy"):
                if (best_metric is None
                        or te_by_PELE_id[PELE_id] < best_metric):
                    best_metric = te_by_PELE_id[PELE_id]
                    best_id = PELE_id

    return ies, rmsds, tes, best_id


def get_ligand_rotatable_bonds(lig_rotamers_path):
    counter = 0
    with open(str(lig_rotamers_path), 'r') as lrl:
        for line in lrl:
            line = line.strip()
            if (line.startswith('sidelib')):
                counter += 1
    return counter


def get_global_metric(metrics_by_PELE_id):
    metrics = []
    for PELE_id, m in metrics_by_PELE_id.items():
        metrics.append(m)

    return np.mean(metrics)


def generate_plot(PELE_ids, filtered_PELE_ids_1, filtered_PELE_ids_2,
                  rmsd_by_PELE_id, energy_by_PELE_id, representative_PELE_id,
                  results, cluster_id,
                  output_path, energy_label):
    #fig = plt.figure()
    ax = plt.subplot(111)

    x = []
    y = []
    for PELE_id in zip(*PELE_ids):
        if (PELE_id in filtered_PELE_ids_1):
            continue
        if (PELE_id[2] == 0):
            continue
        x.append(rmsd_by_PELE_id[PELE_id])
        y.append(energy_by_PELE_id[PELE_id])

    h1 = plt.scatter(x, y, color='grey', alpha=0.5, label='All')

    x = []
    y = []
    for PELE_id in filtered_PELE_ids_1:
        if (PELE_id in filtered_PELE_ids_2):
            continue
        if (PELE_id[2] == 0):
            continue
        x.append(rmsd_by_PELE_id[PELE_id])
        y.append(energy_by_PELE_id[PELE_id])

    h2 = plt.scatter(x, y, color='red', alpha=0.5,
                     label='H bond filter')

    x = []
    y = []
    for r, PELE_id in zip(results, filtered_PELE_ids_2):
        if (r == cluster_id):
            continue
        if (PELE_id[2] == 0):
            continue
        x.append(rmsd_by_PELE_id[PELE_id])
        y.append(energy_by_PELE_id[PELE_id])

    h3 = plt.scatter(x, y, color='green', alpha=0.5,
                     label='Energetic filter')

    x = []
    y = []
    for r, PELE_id in zip(results, filtered_PELE_ids_2):
        if (r == cluster_id):
            if (PELE_id == representative_PELE_id):
                continue
            x.append(rmsd_by_PELE_id[PELE_id])
            y.append(energy_by_PELE_id[PELE_id])

    h4 = plt.scatter(x, y, color='blue', alpha=0.5,
                     label='Selected cluster')

    h5 = plt.scatter(rmsd_by_PELE_id[representative_PELE_id],
                     energy_by_PELE_id[representative_PELE_id], marker='x',
                     color='black', alpha=1, label='Representative\nstructure')

    plt.xlabel('RMSD to initial structure ($\AA$)', fontweight='bold')
    plt.ylabel(energy_label + ' ($kcal/mol$)', fontweight='bold')

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
    plt.legend([h1, h2, h3, h4, h5])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.subplots_adjust(left=0.15, right=0.70)

    plt.savefig(str(output_path.joinpath(
        energy_label.lower().replace(' ', '_') + '_plot.png')))

    plt.close()


def main():
    # Parse args
    PELE_sim_paths, lig_resname, bandwidth, proc_number, \
        ie_col, rmsd_col, topology_relative_path, report_name, \
        hbonds_relative_path, golden_hbonds_1, golden_hbonds_2, \
        minimum_g2_conditions, output_relative_path, generate_plots, \
        representative_extraction_method, maximum_cluster_density_threshold = \
        parse_args()

    golden_hbonds_1 = get_hbond_linkers(golden_hbonds_1)
    golden_hbonds_2 = get_hbond_linkers(golden_hbonds_2)

    check_hbonds_linkers(golden_hbonds_1)
    check_hbonds_linkers(golden_hbonds_2)

    print(' - Golden H bonds set 1:')
    print_hbonds(golden_hbonds_1)
    print(' - Golden H bonds set 2 ({} '.format(minimum_g2_conditions)
          + 'of them need to be fulfilled):')
    print_hbonds(golden_hbonds_2)

    all_sim_it = SimIt(PELE_sim_paths)

    for PELE_sim_path in all_sim_it:
        print('')
        print(' - Extracting ligand coords from {}'.format(PELE_sim_path))
        hbonds_path = PELE_sim_path.joinpath(hbonds_relative_path)
        topology_path = PELE_sim_path.joinpath(topology_relative_path)
        lig_rotamers_path = PELE_sim_path.joinpath('DataLocal/'
                                                   + 'LigandRotamerLibs/'
                                                   + '{}.rot.assign'.format(
                                                       lig_resname))

        if (not topology_path.is_file()):
            print(' - Skipping simulation because topology file with '
                  + 'connectivity was missing')
            continue

        if (not hbonds_path.is_file()):
            print(' - Skipping simulation because hbonds file was '
                  + 'missing')
            continue

        if (not lig_rotamers_path.is_file()):
            print(' - Skipping simulation because ligand rotamer library was'
                  + ' missing')
            continue

        sim_it = SimIt(PELE_sim_path)
        sim_it.build_traj_it('output', 'trajectory', 'xtc')

        trajectories = [traj for traj in sim_it.traj_it]

        print(' - Detected {} trajectories'.format(len(trajectories)))

        reports = get_reports_list(trajectories, report_name)

        print(' - Reading ligand properties')
        ligand_heavy_atoms, ligand_mass = extract_ligand_properties(
            topology_path, lig_resname)
        print('   - Detected {} heavy atoms'.format(ligand_heavy_atoms))
        print('   - Molecular weight: {:.3f}'.format(ligand_mass))

        PELE_ids = extract_PELE_ids(reports)

        hbond_data, _, _ = extract_hbond_linkers(hbonds_path)

        print(' - Detected {} sets of H bonds'.format(len(hbond_data)))

        total_fulfillments, total_models, _, _, f_hbond_data = \
            hbond_fulfillment(hbond_data, golden_hbonds_1, golden_hbonds_2,
                              minimum_g2_conditions)

        if total_fulfillments == 0:
            print(' - Skipping simulation because no H bonds were found')
            continue

        metrics = extract_ligand_metrics(reports, (ie_col, rmsd_col, 4),
                                         proc_number)

        ies = []
        rmsds = []
        tes = []
        for chunk in metrics:
            _ies = []
            _rmsds = []
            _tes = []
            for ie, rmsd, te in chunk:
                _ies.append(float(ie))
                _rmsds.append(float(rmsd))
                _tes.append(float(te))
            ies.append(_ies)
            rmsds.append(_rmsds)
            tes.append(_tes)

        ie_by_PELE_id = get_metric_by_PELE_id(PELE_ids, ies)
        rmsd_by_PELE_id = get_metric_by_PELE_id(PELE_ids, rmsds)
        te_by_PELE_id = get_metric_by_PELE_id(PELE_ids, tes)

        filtered_PELE_ids_1 = filter_by_hbonds(f_hbond_data)

        print(' - H bond filtering: {} '.format(len(filtered_PELE_ids_1))
              + 'structures were selected out of '
              + '{}'.format(len(PELE_ids[0])))

        if (len(filtered_PELE_ids_1) == 0):
            print(' - Skipping simulation because no model fulfills the '
                  + 'required conditions')
            continue

        filtered_PELE_ids_2 = filter_by_energies(filtered_PELE_ids_1,
                                                 ie_by_PELE_id)

        print(' - Energetic filtering: {} '.format(len(filtered_PELE_ids_2))
              + 'structures were selected out of '
              + '{}'.format(len(filtered_PELE_ids_1)))

        if (len(filtered_PELE_ids_2) == 0):
            print(' - Skipping simulation because no model fulfills the '
                  + 'required conditions')
            continue

        lig_coords, filtered_PELE_ids_2 = extract_ligand_coords(
            filtered_PELE_ids_2, trajectories, lig_resname, topology_path,
            proc_number)

        print(' - Clustering filtered ligand coordinates')
        results, cluster_centers = clusterize(lig_coords, bandwidth,
                                              proc_number)

        p_dict = calculate_probabilities(results)

        print(' - Clustering resulted in {} clusters '.format(len(p_dict))
              + 'with densities:')
        for c, f in sorted(p_dict.items(), key=itemgetter(0)):
            print('   - {:3d}: {:3.2f}'.format(c, f))

        cluster_ids, selected_centers = get_most_populated_clusters(
            p_dict, cluster_centers, maximum_cluster_density_threshold)

        cluster_id, ies, rmsds, tes, representative_PELE_id = get_best_cluster(
            cluster_ids, selected_centers, results,
            filtered_PELE_ids_2, ie_by_PELE_id, rmsd_by_PELE_id,
            te_by_PELE_id, lig_coords, representative_extraction_method)

        global_mean_te_energy = get_global_metric(te_by_PELE_id)
        global_mean_ie_energy = get_global_metric(ie_by_PELE_id)

        output_path = PELE_sim_path.joinpath(output_relative_path)

        if (not output_path.is_dir()):
            os.mkdir(str(output_path))

        print(' - Results:')
        print('   - Number of clusters:      {:25d}'.format(
            len(cluster_centers)))
        print('   - Selected cluster:        {:25d}'.format(cluster_id))
        print('   - Mean total energy:       {:25.1f}'.format(np.mean(tes)))
        print('   - Mean interaction energy: {:25.1f}'.format(np.mean(ies)))
        print('   - Mean RMSD (respect to initial structure): '
              + '{:8.1f}'.format(np.mean(rmsds)))
        print('   - Min total energy:        {:25.1f}'.format(np.min(tes)))
        print('   - Min interaction energy:  {:25.1f}'.format(np.min(ies)))
        print('   - Min RMSD (respect to initial structure):  '
              + '{:8.1f}'.format(np.min(rmsds)))
        print('   - Max total energy:        {:25.1f}'.format(np.max(tes)))
        print('   - Max interaction energy:  {:25.1f}'.format(np.max(ies)))
        print('   - Max RMSD (respect to initial structure):  '
              + '{:8.1f}'.format(np.max(rmsds)))
        print('   - Representative structure: '
              + 'epoch: {}, '.format(representative_PELE_id[0])
              + 'trajectory: {}, '.format(representative_PELE_id[1])
              + 'model: {}'.format(representative_PELE_id[2]))

        with open(str(output_path.joinpath('metrics.out')), 'w') as f:
            f.write('{:19}: {:10d}\n'.format(
                'Heavy atoms', ligand_heavy_atoms))
            f.write('{:19}: {:10.3f}\n'.format(
                'Molecular weight', ligand_mass))
            f.write('{:19}: {:10d}\n'.format(
                'N. rotatable bonds',
                get_ligand_rotatable_bonds(lig_rotamers_path)))
            f.write('{:19}: {:10d}\n'.format(
                'N. clusters', len(cluster_centers)))
            f.write('{:19}: {: 10.1f}\n'.format(
                'Global Mean Tot. E.', global_mean_te_energy))
            f.write('{:19}: {: 10.1f}\n'.format(
                'Global Mean Int. E.', global_mean_ie_energy))

        with open(str(output_path.joinpath('cluster.out')), 'w') as f:
            f.write('{:12s}: {: 10d}\n'.format('Cluster size', len(tes)))
            f.write('{:12s}: {: 10.1f}\n'.format('Mean Tot. E.', np.mean(tes)))
            f.write('{:12s}: {: 10.1f}\n'.format('Mean Int. E.', np.mean(ies)))
            f.write('{:12s}: {: 10.1f}\n'.format('Mean RMSD', np.mean(rmsds)))
            f.write('{:12s}: {: 10.1f}\n'.format('Min Tot. E.', np.min(tes)))
            f.write('{:12s}: {: 10.1f}\n'.format('Min Int. E.', np.min(ies)))
            f.write('{:12s}: {: 10.1f}\n'.format('Min RMSD', np.min(rmsds)))
            f.write('{:12s}: {: 10.1f}\n'.format('Max Tot. E.', np.max(tes)))
            f.write('{:12s}: {: 10.1f}\n'.format('Max Int. E.', np.max(ies)))
            f.write('{:12s}: {: 10.1f}\n'.format('Max RMSD', np.max(rmsds)))

        rep_traj = md.load(str(PELE_sim_path.joinpath(
            'output/{}/trajectory_{}.xtc'.format(
                *representative_PELE_id[:2]))),
            top=str(topology_path))

        rep_model = rep_traj[representative_PELE_id[2]]
        rep_model.save_pdb(str(output_path.joinpath('rep_structure.pdb')))

        if (generate_plots):
            generate_plot(PELE_ids, filtered_PELE_ids_1, filtered_PELE_ids_2,
                          rmsd_by_PELE_id, ie_by_PELE_id,
                          representative_PELE_id, results, cluster_id,
                          output_path, 'Interaction energy')
            generate_plot(PELE_ids, filtered_PELE_ids_1, filtered_PELE_ids_2,
                          rmsd_by_PELE_id, te_by_PELE_id,
                          representative_PELE_id, results, cluster_id,
                          output_path, 'Total energy')


if __name__ == "__main__":
    main()
