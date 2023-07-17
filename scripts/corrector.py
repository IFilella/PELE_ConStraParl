# -*- coding: utf-8 -*-
"""
This module is designed to correct energies and obtain
strain data of a simulation.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
import shutil

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from seaborn import displot
from collections import OrderedDict
from collections import defaultdict

# ----------------------------------------------------------------------- #
# Constant:
T = 298.
R = 1.985e-3
# ----------------------------------------------------------------------- #


def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.
    Parameters

    Parameters
    ----------
    - args : list[str]
        List of command-line arguments to parse

    Returns
    ----------
    - parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='LIG_Pele', help="Name of the directory where the simulation\
        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files used for the simulation.")
    parser.add_argument("-q", "--quantile", type=float, dest="quantile",
                        default=1, help="Percentage of data with lowest total energy snapshots you want\
        to keep to assess the strain energy")
    parser.add_argument("-cl", "--clusters_folder", type=str, dest="clusters_folder",
                        default=None, help="Name of the folder you want to obtain the clusters information from.")
    parser.add_argument("-lf", "--ligand_folder", type=str, dest="ligand_folder",
                        default='LIG_ligand', help="Name of the output folder of the ligand simulation.")

    parser.add_argument("--skip_strain_per_cluster", dest="strain_per_cluster_bool",
                        default=False, action='store_true', help="Flag to choose if strain per cluster code is skipped.")
    parser.add_argument("--skip_platform_analysis", dest="analysis_bool",
                        default=False, action='store_true', help="Flag to choose if the platform analysis is skipped.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def corrector(input_folder,
              ligand_folder,
              residue_name,
              report_name,
              quantile,
              clusters_folder,
              strain_per_cluster_bool,
              analysis_bool):
    """
    Function
    ----------
    With all the corrections calculated, they are introduced in the reports.

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.
    - residue_name : str
        Residue name of the ligand in the pdb of each cluster
    -report_name : str
        Name of the report files we want to correct.
    - quantile : float
        Percentage (0,1] of data with lowest total energy snapshots you want
        to keep to assess the strain energy
    """

    def path_definer(input_folder,
                     ligand_folder,
                     residue_name):
        """
        Function
        ----------
        Defines all the paths that are going to be used

        Parameters
        ----------
        - input_folder : str
            The path to the directory created by the induced fit simulation.
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster

        Returns
        ----------
        - path_pl_simulation: str
            The path to the protein-ligand simulation.
        - path_pl_output : str
            The path to the protein-ligand simulation output.
        - path_l_simulation : str
            The path to the ligand simulation.
        - path_pl_results : str
            The path to the folder where results will be stored.
        """

        path = str(pathlib.Path().absolute())

        path_pl_simulation = os.path.join(path, input_folder)
        path_pl_output = os.path.join(path_pl_simulation, 'output')
        path_pl_results = os.path.join(path_pl_simulation, 'strain')

        path_l_simulation = os.path.join(path, ligand_folder)

        if os.path.isdir(path_pl_simulation) == False:
            raise Exception('PathError: There is no folder with this path: ' +
                            path_pl_simulation + '. Please check the path and the folder name.')

        if os.path.isdir(path_pl_results) == False:
            os.mkdir(path_pl_results)

        return path_pl_simulation, path_pl_output, path_l_simulation, path_pl_results

    def corrections_detector(path_pl_simulation,
                             path_l_simulation):
        """
        Function
        ----------
        Checks whether there is a correction to be implemented to the reports of the
        protein ligand simulation.

        Parameters
        ----------
        - path_pl_simulation : str
            The path to the protein-ligand simulation.
        - path_l_simulation : str
            The path to the ligand simulation.

        Returns
        ----------
        - correction_number : int
            Number that indicates univocally what corrections are to be implemented.

        """

        correction_number = 0

        if os.path.isfile(os.path.join(path_l_simulation, 'energy.csv')) == True:

            print(' -   Strain correction found.')
            correction_number += 1

        else:

            print(' -   Strain correction couldn\'t be found.'
                  '     File containing energy information (energy.csv)\n'
                  '     has not been found in\n\n'
                  '    ' + path_l_simulation + '.\n')

        if os.path.isfile(os.path.join(path_pl_simulation, 'entropy.csv')) == True:

            print(' -   Entropy correction found.')
            correction_number += 2

        else:

            print(' -   Entropy correction couldn\'t be found.\n'
                  '     File containing entropy information (entropy.csv)\n'
                  '     has not been found in\n')
            print('     ' + path_pl_simulation)
            print(' ')

        if correction_number == 0:

            raise Exception('NoCorrectionToImplement: Either there is no correction to apply to the data' +
                            ' or they have not been found in the paths they are supposed to. (1) ' + path_pl_simulation + ' (2) ' + path_l_simulation + '.')

        return correction_number

    def ligand_min(path_l_simulation):
        """
        Function
        ----------
        Retrieves the energetic scoring of the ligand (we have chosen).

        Parameters
        ----------
        - path_l_simulation : str
            The path to the ligand simulation.

        Returns
        ----------
        - ligand_min_energy : float
            Ligand's minimum energy of the scoring function chosen.
        """

        cont = 0

        with open(os.path.join(path_l_simulation, 'energy.csv')) as filein:

            for line in filein:

                if cont != 0:

                    # This could change depending on the user.
                    line = line.split(',')
                    ligand_min_energy = float(line[0])

                cont += 1

        return ligand_min_energy

    def entropy_correction(path_pl_simulation):
        """
        Function
        ----------
        Retrieves the entropic change of the ligand.

        Parameters
        ----------
        - path_pl_simulation : str
            The path to the protein-ligand simulation.

        Returns
        ----------
        - entropy_change : float
            Entropy change calculated from the protein-ligand and the ligand simulations.
        """

        cont = 0

        with open(os.path.join(path_pl_simulation, 'entropy.csv')) as filein:

            for line in filein:

                if cont != 0:

                    line = line.split(',')
                    entropy_change = float(line[3])

                cont += 1

        return entropy_change

    def column_retriever(file):
        """
        Function
        ----------
        Retrieves the position of the binding energy.

        Parameters
        ----------
        - file : list
            The path to the first report.

        Returns
        ----------
        - column_current_energy : int
            Column where the interesting total energy is located in the report.
        - column_binding_energy : int
            Column where the interesting binding energy is located in the report.
        - column_internal_energy : int
            Column where the interesting internal energy is located in the report.
        """

        cont = 0

        with open(file, 'r') as filein:

            for line in filein:

                if cont == 0:

                    line = line.split('    ')
                    column_current_energy = line.index('currentEnergy') + 1
                    column_binding_energy = line.index('BindingEnergy') + 1
                    column_internal_energy = line.index('InternalEnergy') + 1

                    #
                    print(
                        ' -   Column for the binding energy:', column_binding_energy - 1)
                    print(
                        ' -   Column for the internal energy:', column_internal_energy - 1)
                    #

                cont += 1

        return column_current_energy, column_binding_energy, column_internal_energy

    def trajectory_retriever(path):
        """
        Function
        ----------
        Retrieves the file format of trajectories.

        Parameters
        ----------
        - path : str
            The path to the epoch.

        Returns
        ----------
        - file_format : int
            Format of the trajectories.
        """

        trajectory_list = [x for x in os.listdir(path) if x.startswith('trajectory_')]
        file_format = trajectory_list[0].split('trajectory_')[-1].split('.')[-1]

        return file_format

    def correction_implementer(column_current_energy,
                               column_binding_energy,
                               column_internal_energy,
                               path_pl_output,
                               report_name,
                               ligand_min_energy,
                               entropy_change):
        """
        Function
        ----------
        Correct the reports with the calculated corrections at disposal and
        store simulation data.

        Parameters
        ----------
        - column_current_energy : int
            Column where the interesting total energy is located in the report.
        - column_binding_energy : int
            Column where the interesting binding energy is located in the report.
        - column_internal_energy : int
            Column where the interesting internal energy is located in the report.
        - path_pl_output : str
            The path to the protein-ligand simulation output.
        - report_name : str
            Name of the report files we want to correct.
        - ligand_min_energy : float
            Ligand's minimum energy of the scoring function chosen.
        - entropy_change : float
            Entropy change calculated from the protein-ligand and the ligand simulations.

        Returns
        ----------
        - strain_energy_list : list
            List of all the strain energies calculated.
        - simulation_df : pd.DataFrame
            Data frame with currentEnergy and strainEnergy values of the simulation.
        """

        def copier(path_pl_output,
                   report_name):
            """
            Function
            ----------
            Copies al the reports in place with the prefix mod_ added. It also
            stores information for the correction.

            Parameters
            ----------
            - path_pl_output : str
                The path to the protein-ligand simulation output.
            - report_name : str
                Name of the report files we want to correct.

            Returns
            ----------
            - report_paths : list
                List of all the reports' paths.
            """

            cont_reports = 0
            report_paths = []

            # Copying reports
            if os.path.isdir(path_pl_output):

                files = os.listdir(path_pl_output)

                for folder in files:

                    if folder.isnumeric():

                        full_path = os.path.join(path_pl_output, folder)

                        files_subdir = os.listdir(full_path)

                        for report in files_subdir:

                            if report.startswith(report_name) and report.split(report_name + '_')[1].isnumeric():

                                shutil.copy(os.path.join(full_path, report),
                                            os.path.join(full_path, 'mod_' + report_name + '_' + report.split(report_name + '_')[1]))

                                report_paths.append(
                                    os.path.join(full_path, report))
                                cont_reports += 1

                            else:
                                continue

                        if cont_reports == 0:

                            raise Exception('ReportNameError: No reports beginning with \"' + report_name + '\" were found in '
                                            + full_path + '.')

                    else:

                        continue

            return report_paths

        def corrector(column_binding_energy,
                      column_internal_energy,
                      report_paths,
                      report_name,
                      ligand_min_energy,
                      entropy_change):
            """
            Function
            ----------
            Correct the new reports with the calculated corrections.

            Parameters
            ----------
            - column_binding_energy : int
                Column of the report where the binding energy metric is located.
            - column_internal_energy : int
                Column of the report where the internal energy metric is located.
            - report_paths : list
                List of all the reports' paths.
            - report_name : str
                Name of the report files we want to correct.
            - ligand_min_energy : float
                Ligand's minimum energy of the scoring function chosen.
            - entropy_change : float
                Entropy change calculated from the protein-ligand and the ligand simulations.

            Returns
            ----------
            - strain_energy_list : list
                List with all the strain energies calculated in a simulation.
            - simulation_df : pd.DataFrame
                Data frame with current energies and strain information.
            """

            epochs = []
            trajectories = []
            models = []
            currentEnergies = []
            strainEnergies = []

            for report_path in report_paths:

                epoch = report_path.split('/')[-2]
                trajectory = report_path.split('/')[-1].split('_')[-1]

                path_modified_report = report_path.replace(
                    report_name, 'mod_' + report_name)

                with open(path_modified_report, 'w') as fileout:

                    cont = 0

                    with open(report_path) as filein:

                        for line in filein:

                            if cont == 0:

                                line = line.split()
                                line.append('strainEnergy')
                                fileout.write(
                                            "    ".join(str(v) for v in line) + '\n')

                            else:

                                if '/0/' + report_name in report_path:

                                    if cont == 1:

                                        line = line.split()
                                        strain_energy = float(
                                            line[column_internal_energy-1]) - ligand_min_energy
                                        line[column_binding_energy-1] = \
                                            str("{:.4f}".format(float(
                                                line[column_binding_energy-1]) + strain_energy + entropy_change))
                                        line.append("{:.4f}".format(strain_energy))
                                        fileout.write(
                                            "     ".join(str(v) for v in line) + '\n')

                                    else:

                                        line = line.split()

                                        currentEnergy = line[column_current_energy-1]
                                        model = line[2]
                                        strain_energy = float(
                                            line[column_internal_energy-1]) - ligand_min_energy

                                        line[column_binding_energy-1] = \
                                            str("{:.4f}".format(float(
                                                line[column_binding_energy-1]) + strain_energy + entropy_change))
                                        line.append("{:.4f}".format(strain_energy))
                                        fileout.write(
                                            "     ".join(str(v) for v in line) + '\n')

                                        epochs.append(int(epoch))
                                        trajectories.append(int(trajectory))
                                        models.append(int(model))
                                        currentEnergies.append(
                                            float(currentEnergy))
                                        strainEnergies.append(
                                            float(strain_energy))

                                else:

                                    line = line.split()

                                    currentEnergy = line[column_current_energy-1]
                                    model = line[2]
                                    strain_energy = float(
                                        line[column_internal_energy-1]) - ligand_min_energy

                                    line[column_binding_energy-1] = \
                                        str("{:.4f}".format(float(
                                            line[column_binding_energy-1]) + strain_energy + entropy_change))
                                    line.append("{:.4f}".format(strain_energy))
                                    fileout.write("     ".join(str(v) for v in line) + '\n')

                                    epochs.append(int(epoch))
                                    trajectories.append(int(trajectory))
                                    models.append(int(model))
                                    currentEnergies.append(
                                        float(currentEnergy))
                                    strainEnergies.append(float(strain_energy))

                            cont += 1

            simulation_df = pd.DataFrame(OrderedDict([('epoch', epochs),
                                                      ('trajectory', trajectories),
                                                      ('model', models),
                                                      ('currentEnergy',
                                                       currentEnergies),
                                                      ('strainEnergy', strainEnergies)]))

            simulation_df = simulation_df.sort_values(
                ['epoch', 'trajectory', 'model', 'currentEnergy', 'strainEnergy'], ascending=True).reset_index(drop=True)

            return strainEnergies, simulation_df

        # Copying  reports
        report_paths = copier(path_pl_output, report_name)

        # Correcting copied reports and storing information
        strain_energy_list, simulation_df = corrector(column_binding_energy,
                                                      column_internal_energy,
                                                      report_paths,
                                                      report_name,
                                                      ligand_min_energy,
                                                      entropy_change)

        return strain_energy_list, simulation_df

    def quantile_retriever(simulation_df,
                           quantile):
        """
        Function
        ----------
        Obtain a list of strain energies corresponding to those conformations
        with energies below a certain quantile.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Data frame with current energies and strain information.
        - quantile : float
            Percentage of data with lowest total energy snapshots you want to keep
            to assess the strain energy.

        Returns
        ----------
        - strain_energy_quantile : list
            List with the strain energies belonging to the quantile percentage of data
            with lowest total energy.
        """

        quantile_value = simulation_df.currentEnergy.quantile(quantile)
        simulation_df_quantile = simulation_df[simulation_df.currentEnergy < quantile_value]
        strain_energy_quantile = list(
            np.array(simulation_df_quantile['strainEnergy']))

        return strain_energy_quantile

    def strain_per_cluster(clusters_folder,
                           simulation_df,
                           path_pl_simulation,
                           path_pl_results):
        """
        Function
        ----------
        Obtaining average strain per cluster obtained from the platform.

        Parameters
        ----------
        - clusters_folder : str
            Name of the clusters folder you want to obtain the information
            from.
        - simulation_df : pd.DataFrame
            Data frame with current energies and strain information.
        - path_pl_simulation : str
            The path to the protein-ligand simulation.
        - path_pl_results : str
            The path to the folder where results will be stored.
        """

        def cluster_information_retriever(path_cluster,
                                          path_pl_results):
            """
            Function
            ----------
            Calculating average strain per cluster with simulation_df strain
            information and data.csv cluster information from the directory
            generated by the platform.

            Parameters
            - path_cluster : str
                The path to the location selected to obtain the cluster
                information.
            - path_pl_results : str
                The path to the folder where results will be stored.
            """

            results = defaultdict(list)

            with open(os.path.join(path_cluster, 'data.csv')) as filein:

                cont = 0

                for line in filein:

                    if cont != 0:

                        line = line.split(',')
                        cluster = line[-1]

                        try:

                            cluster = int(cluster)
                            epoch = int(line[-3])
                            trajectory = int(
                                line[-2].split('/')[-1].split('_')[1].split('.pdb')[0])
                            model = int(line[1])

                            print(simulation_df)

                            epoch_df = simulation_df.loc[simulation_df['epoch'] == epoch]
                            epoch_trajectory_df = epoch_df.loc[epoch_df['trajectory']
                                                               == trajectory]
                            epoch_trajectory_model_df = epoch_trajectory_df.loc[
                                epoch_trajectory_df['model'] == model]

                            value = epoch_trajectory_model_df['strainEnergy'].values[0]
                            results[cluster].append(value)

                        except ValueError:
                            pass

                    cont += 1

            cluster_strain = list(results.items())
            cluster_strain_dict = {}

            for cluster_num, strains in cluster_strain:

                average_strain = np.average(np.array(strains))
                cluster_strain_dict[cluster_num] = [average_strain]

            cluster_strain_df = pd.DataFrame.from_dict(cluster_strain_dict)
            cluster_strain_df.rename(
                index={0: 'strain (kcal/mol)'}, inplace=True)
            cluster_strain_df.index.name = 'clusters'
            cluster_strain_df.to_csv(os.path.join(
                path_pl_results, 'strain_cluster.csv'))

        if clusters_folder is None:

            path_cluster_analysis = os.path.join(
                path_pl_simulation, 'analysis')
            path_cluster_results = os.path.join(path_pl_simulation, 'results')

            if os.path.isdir(path_cluster_analysis) == False:

                if os.path.isdir(path_cluster_results) == False:

                    #
                    print('     -   No /analysis or /results folder has been found in\n\n   ' +
                          path_pl_simulation + '\n      Strain per cluster cannot be calculated.')
                    #

                else:

                    #
                    print('     -   Information will be obtained from /results.')
                    #

                    cluster_information_retriever(path_cluster_results,
                                                  path_pl_results)

                    #
                    print('     -   Strain per cluster stored in results folder.\n')
                    #

            else:

                #
                print('     -   Information will be obtained from /analysis.')
                #

                cluster_information_retriever(path_cluster_results,
                                              path_pl_results)

                #
                print('     -   Strain per cluster stored in ' +
                      input_folder + '/strain/strain_cluster.csv.\n')

        else:

            #
            print('     -   Information will be obtained from /' +
                  clusters_folder + '.\n')
            #

            path_cluster = os.path.join(path_pl_simulation, clusters_folder)

            cluster_information_retriever(path_cluster,
                                          path_pl_results)

    def boltzmann_weighted(be, ene_t, T):
        """
        Function
        ----------
        Calculates boltzmann weighted energy.

        Parameters
        ----------
        - be : list
            Binding energies of all the simulation.
        - ene_t : list
            Total energies of all the simulation.
        - T : float
            Temperature of the simulation.

        Returns
        ----------
        - ene_bz : float
            Value of the boltzmann weighted energy.
        """

        exp_bz = np.exp(-ene_t/(R*T))
        nominator = be.dot(exp_bz)
        denominator = np.sum(exp_bz)
        ene_bz = nominator/denominator

        return ene_bz

    def analysis_files_writer(column_binding_energy,
                              file_format,
                              path_pl_simulation):
        """
        Function
        ----------
        Write files for further analysis with pele platform.

        Parameters
        ----------
        - column_binding_energy : int
            Column of the report where the binding energy metric is located.
        - file_format : str
            Format of the trajectory files.
        - path_pl_simulation: str
            The path to the protein-ligand simulation.

        Returns
        ----------
        - path_to_run : str
            Path to the slurm file to be run.
        """

        path_to_run = os.path.join(path_pl_simulation, 'run_analysis')
        list_input_files = [x for x in os.listdir(os.path.join(path_pl_simulation,'input'))]

        if 'receptor.pdb' in list_input_files:
            list_input_files.remove('receptor.pdb')
        if 'ligand.pdb' in list_input_files:
            list_input_files.remove('ligand.pdb')

        #topology_file = os.path.join(path_pl_simulation,'input',list_input_files[0])
	topology_file = os.path.join(path_pl_simulation,'output','topologies','topologies_0.pdb')
        with open(path_to_run, 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=analysis\n'
                '#SBATCH --output=analysis.out\n'
                '#SBATCH --error=analysis.err\n'
                '#SBATCH --ntasks=1\n'
                '#SBATCH --qos=bsc_ls\n'
                '#SBATCH --time=00-00:30:00\n'
                '\n'
                'module load ANACONDA/2019.10\n'
                'module load intel mkl impi gcc # 2> /dev/null\n'
                'module load impi\n'
                'module load boost/1.64.0\n'
                '\n'
                'eval "$(conda shell.bash hook)"\n'
                'conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.3\n'
                '\n'
                'python script.py\n'
            )

        with open(os.path.join(path_pl_simulation, 'script.py'), 'w') as fileout:

            fileout.writelines(
                'from pele_platform.analysis import Analysis\n'
                '\n'
                'analysis = Analysis(resname="' + residue_name +
                '", chain="L", simulation_output="output", be_column = ' + str(column_binding_energy) + ', report="' +
                'mod_' + report_name + '", traj="trajectory.' + file_format + '",  topology = "' + topology_file + '", cpus=48)\n'
                'analysis.generate(path="analysis", clustering_type="meanshift")\n'
            )

        return path_to_run

    def results_writer(strain_energy_list,
                       strain_energy_quantile,
                       path,
                       report_name,
                       quantile,
                       ligand_min_energy,
                       strain_bz):
        """
        Function
        ----------
        Writes and plots result files of strain values.

        Parameters
        ----------
        - strain_energy_list : list
            List with all the strain values calculated from the entire simulation.
        - strain_energy_quantile : list
            List with strain values of the entire simulation corresponding to a certain quantile.
        - path : str
            The path to the protein-ligand strain results folder.
        - report_name : str
            Name of the report files we want to correct.
        - quantile : float
            Percentage of data with lowest total energy snapshots you want to keep
            to assess the strain energy.
        - ligand_min_energy : float
            Ligand's minimum energy of the scoring function chosen.
        - strain_bz : float
            Boltzmann averaged strain energy of the simulation.
        """

        def histogram_function(strain_energy):
            """
            Function
            ----------
            Writes and plots result files of strain values.

            Parameters
            ----------
            - strain_energy : list
                List with all the strain values calculated from the entire simulation.
            """

            strain_energy = [
                item for item in strain_energy if item < 200]

            strain_energy_vector = np.array(strain_energy)

            bin_edges = np.histogram_bin_edges(strain_energy_vector, bins='fd')
            density, _ = np.histogram(strain_energy_vector, bins=bin_edges)

            hist_ene = 0.5*(bin_edges[np.argmax(density)] +
                            bin_edges[np.argmax(density) + 1])
            minimum_ene = min(strain_energy_vector)
            average_ene = np.average(strain_energy_vector)
            max_ene = max(strain_energy_vector)

            return strain_energy_vector, bin_edges, hist_ene, minimum_ene, average_ene, max_ene

        if quantile == 1.:

            strain_energy_vector, bin_edges, hist_ene, minimum_ene, average_ene, max_ene = histogram_function(
                strain_energy_list)

            if min(strain_energy_vector) < 0.:

                strain_energy_vector = strain_energy_vector - minimum_ene

                #
                print('\n'
                      '                              WARNING:                               \n'
                      '   Lower ligand energies were found in the induced fit simulations.  \n'
                      '   The results in mod_' + report_name + ' have been corrected with   \n'
                      '   Ligand minimum = ' + "{:3.3f}".format(ligand_min_energy) + ' > ' +
                      "{:3.3f}".format(minimum_ene + ligand_min_energy) +
                      ' = Induced fit minimum.\n'
                      '   Better sampling of the ligand is required.          \n'
                      )
                #

            with open(os.path.join(path, 'strain.csv'), 'w') as fileout:
                fileout.writelines(
                    'Minimum,Histogram max,Average,Maximum,Boltzmann\n'
                    '' + str(minimum_ene) + ',' + str(hist_ene) + ',' +
                    str(average_ene) + ',' + str(max_ene) + ',' + str(strain_bz) + '\n')

            # Plot
            displot(strain_energy_vector, kind="kde", color='black', label='KDE plot')
            plt.title('Strain distribution')
            plt.hist(strain_energy_vector, bins=bin_edges,
                     density=True, color='#B3C6CF')
            plt.axvline(x=hist_ene, color = '#0B3954', label='Hist max = ' + str("{:.3f}".format(hist_ene)))
            plt.axvline(x=average_ene, color = '#CE796B', label='Average = ' + str("{:.3f}".format(average_ene)))
            plt.legend(loc="best")
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.tight_layout()
            plt.savefig(os.path.join(path, 'density_strain.svg'), format='svg', bbox_inches = "tight", transparent=True)

        else:

            print(' -   Quantile chosen for the results:', quantile, '\n')

            strain_energy_vector, \
                bin_edges, \
                hist_ene,\
                minimum_ene,\
                average_ene,\
                max_ene = histogram_function(strain_energy_list)

            strain_energy_vector_q,\
                bin_edges_q,\
                hist_ene_q,\
                minimum_ene_q,\
                average_ene_q,\
                max_ene_q = histogram_function(strain_energy_quantile)

            with open(os.path.join(path, 'strain' + str(int(quantile*100)) + '.csv'), 'w') as fileout:
                fileout.writelines(
                    'Minimum,Histogram max,Average,Maximum\n'
                    '' + str(minimum_ene) + ',' + str(hist_ene) + ',' +
                    str(average_ene) + ',' + str(max_ene) +
                    ',' + str(strain_bz) + '\n'
                    '' + str(minimum_ene_q) + ',' + str(hist_ene_q) + ',' +
                    str(average_ene_q) + ',' + str(max_ene_q) + '\n'
                )

            # Plot
            plt.title('Strain distribution')
            plt.hist(strain_energy_vector, bins=bin_edges,
                     density=True, color="blue", label='original')
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.savefig(os.path.join(path, 'density_strain.png'), format='png')

            plt.title('Strain distribution quantile ' +
                      str(int(quantile*100)) + '%')
            plt.hist(strain_energy_vector_q, bins=bin_edges_q,
                     density=True, color="red", label='quantile')
            plt.xlabel('Strain (kcal/mol)')
            plt.ylabel('Density')
            plt.legend(loc='best')
            plt.savefig(os.path.join(path, 'density_strain_quantile' +
                        str(int(quantile*100)) + '.png'), format='png')

    #
    print(' ')
    print('*******************************************************************')
    print('*                        peleCorrector                            *')
    print('* --------------------------------------------------------------- *')
    print('*                   Corrector of induced fit results              *')
    print('*******************************************************************')
    print(' ')
    #

    path_pl_simulation,\
        path_pl_output,\
        path_l_simulation,\
        path_pl_results = path_definer(input_folder, ligand_folder, residue_name)

    correction_number = corrections_detector(
        path_pl_simulation, path_l_simulation)

    # Corrections values ---
    if correction_number == 1:

        ligand_min_energy = ligand_min(path_l_simulation)
        entropy_change = 0.

    elif correction_number == 2:

        ligand_min_energy = 0.
        entropy_change = entropy_correction(path_pl_simulation)

    elif correction_number == 3:

        ligand_min_energy = ligand_min(path_l_simulation)
        entropy_change = entropy_correction(path_pl_simulation)

    # Corrections column location ---
    file = os.path.join(path_pl_simulation, 'output', '1', report_name + '_1')

    column_current_energy, column_binding_energy, column_internal_energy = \
        column_retriever(file)

    # Trajectories format
    file_format = trajectory_retriever(os.path.join(path_pl_simulation, 'output', '1'))

    #
    print(' ')
    print(' -   Implementing corrections...')
    #

    strain_energy_list, simulation_df = correction_implementer(column_current_energy,
                                                               column_binding_energy,
                                                               column_internal_energy,
                                                               path_pl_output,
                                                               report_name,
                                                               ligand_min_energy,
                                                               entropy_change)

    strain_energy_quantile = quantile_retriever(simulation_df,
                                                quantile)

    #
    print('     -   Job finished succesfully. Energies corrected will be found in \n'
          '         mod_' + report_name + ' files.\n')
    #

    if not strain_per_cluster_bool:

        #
        print(' -   Associating strain to RMSD clustered positions.')
        #

        strain_per_cluster(clusters_folder,
                           simulation_df,
                           path_pl_simulation,
                           path_pl_results)

    else:

        #
        print(' -   Strain association to clustered positions skipped.\n')
        #

    strain_bz = boltzmann_weighted(simulation_df['strainEnergy'].to_numpy(),
                                   simulation_df['currentEnergy'].to_numpy(
    ) - min(simulation_df['currentEnergy'].tolist()),
        T)

    results_writer(strain_energy_list,
                   strain_energy_quantile,
                   path_pl_results,
                   report_name,
                   quantile,
                   ligand_min_energy,
                   strain_bz)

    path_to_run = analysis_files_writer(column_binding_energy,
                                        file_format,
                                        path_pl_simulation)

    #
    print(' -   run_analysis and script.py files have been generated.')
    print(' -   All results are stored in ../' + input_folder + '/strain.\n')
    #

    if not analysis_bool:

        #
        print(' -   Launching platform job:\n')
        #

        os.chdir(path_pl_simulation)
        os.system("sbatch %s" % path_to_run)

        #
        print(' ')

    else:

        #
        print(' -   Platform analysis skipped.')
        print(' -   It can be run afterwards with: sbatch run_analysis \n')


def main(args):
    """
    Function
    ----------
    It reads the command-line arguments and runs lice_results.

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    corrector(input_folder=args.input_folder,
              ligand_folder=args.ligand_folder,
              residue_name=args.residue_name,
              report_name=args.report_name,
              quantile=args.quantile,
              clusters_folder=args.clusters_folder,
              strain_per_cluster_bool=args.strain_per_cluster_bool,
              analysis_bool=args.analysis_bool)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)


