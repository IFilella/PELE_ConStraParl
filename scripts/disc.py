# -*- coding: utf-8 -*-
"""
This module is designed to calculate different scoring methods
for the binding energy of a simulation. Different Scorings Calculator
(DiSC).
"""

# Imports
import sys
import os
import pathlib
import argparse
import numpy as np
import time
import matplotlib.pyplot as plt


# ----------------------------------------------------------------------- #
# Constants:
T = 298.
R = 1.985e-3
# ----------------------------------------------------------------------- #


def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.

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
                        default='output', help="Name of the output directory where the files to\
        analyze are located.")
    parser.add_argument("-rn", "--report_name", type=str, dest="report_name",
                        default='report', help="Name of the report files.")
    parser.add_argument("-c", "--column", type=int, dest="column",
                        default=None, help="Column of the metric to consider.")
    parser.add_argument("-T", "--temperature", type=float, dest="temperature",
                        default=298., help="Temperature of the experiment.")
    parser.add_argument("-a", "--action", type=str, dest="action",
                        default='all', help="Function the user wants the script to do: all or evolution.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def statistics(input_folder,
               report_name,
               column,
               T,
               action):
    """
    Function
    ----------
    Reads the data from all the reports and calculates different scoring functions.
    The information is stored in energy.csv

    Parameters
    ----------
    - input_folder : str
        The path to the directory where the output of the PELE simulation is located.
    - report_name : str
        Name of the reports containing the energetic data of all the simulation.
    - column : int
        Column where the interesting metric is located inside the report.
    - T : float
        Temperature to perform the Boltzmann weights with.
    - action : str
        Functionality of the script we want to use.
    """

    def path_definer(input_folder):
        """
        Function
        ----------
        Defines the important paths that are going to be used throughout the simulation.

        Parameters
        ----------
        - input_folder : str
            Name of the folder where the output of the simulation is located.

        Returns
        ----------
        - path_output : str
            Path to the output folder of the simulation.
        """

        path = str(pathlib.Path().absolute())
        path_output = os.path.join(path, input_folder)

        return path_output, path

    def reader(folderpath,
               report_name):
        """
        Function
        ----------
        Reads the data from all the reports and stores it in lists.

        Parameters
        ----------
        - files : list
            The path to the directory where the output of the PELE simulation is located.
        - folderpath : str
            Path where the different epochs of the simulation are located.
        - report_name : str
            Name of the reports to obtain the data from


        Returns
        ----------
        - be : list
            Binding energies of all the simulation.
        - te : list
            Total energies of all the simulation.
        - step : list
            Steps associated to poses for all the simulation.
        - column_be : int
            Column where the Binding energy is located.
        """

        def file_reader(files,
                        folderpath,
                        report_name,
                        column_be,
                        column_te):
            """
            Function
            ----------
            Reads the data from all the reports in a list.

            Parameters
            ----------
            - files : list
                The path to the directory where the output of the PELE simulation is located.
            - folderpath : str
                Path where the different epochs of the simulation are located.
            - report_name : str
                Name of the reports to obtain the data from
            - column_be : int
                Column where the binding energy data is located.
            - column_be : int
                Column where the total energy data is located.

            Returns
            ----------
            - be : list
                Binding energies of all the simulation.
            - te : list
                Total energies of all the simulation.
            - step : list
                Steps associated to poses for all the simulation.
            """

            new_directory = folderpath

            for file in files:

                cont = 0

                if file.startswith(report_name):

                    file_path = os.path.join(new_directory, file)

                    with open(file_path, 'r') as filein:

                        for line in filein:

                            line = line.split('    ')

                            if cont == 1:

                                be.append(float(line[column_be-1]))
                                te.append(float(line[column_te-1]))
                                step.append(int(line[1]))

                            elif cont > 1:

                                repeat = int(line[1]) - step[-1]

                                step.append(int(line[1]))
                                be.extend([float(line[column_be-1])
                                          for _ in range(repeat)])
                                te.extend([float(line[column_te-1])
                                          for _ in range(repeat)])

                            cont += 1

            return be, te, step

        def column_retriever(file):
            """
            Function
            ----------
            Retrieves the position of the binding energy and current energy.

            Parameters
            ----------
            - file : list
                The path to the one report.

            Returns
            ----------
            - column_be : int
                Column where the binding energy data is located.
            - column_te : int
                Column where the total energy data is located.
            """

            cont = 0

            with open(file, 'r') as filein:

                for line in filein:

                    if cont == 0:

                        line = line.split('    ')

                        if 'bindingEnergy' in line:

                            column_be = line.index('bindingEnergy') + 1

                        else:

                            column_be = 1

                            #
                            print(
                                ' -   No bindingEnergy column.')
                            #

                        column_te = line.index('currentEnergy') + 1

                    cont += 1

            return column_be, column_te

        be = []
        te = []
        step = []
        files = os.listdir(folderpath)
        numeric_files = [s for s in os.listdir(folderpath) if s.isnumeric()]

        if len(numeric_files) != 0:

            column_file = os.path.join(
                folderpath, numeric_files[0], report_name + '_1')
            column_be, column_te = column_retriever(column_file)

            for document in numeric_files:

                new_directory = os.path.join(folderpath, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)

                    if report_name in files == False:
                        raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                        Please check the path to the files and the files name.')

                    if column is None:
                        be, te, step = file_reader(files,
                                                   new_directory,
                                                   report_name,
                                                   column_be,
                                                   column_te)

                    else:
                        be, te, step = file_reader(files,
                                                   new_directory,
                                                   report_name,
                                                   column,
                                                   column_te)


        else:

            if report_name in files == False:
                raise Exception('FilePathError: There is no file containing ' + report_name + ' in it. \
                Please check the path to the files and the files name.')

            column_file = os.path.join(folderpath, report_name + '_1')
            column_be, column_te = column_retriever(column_file)

            if column is None:

                if column_be == 1:

                    #
                    print(
                        ' -   Picking currentEnergy as metric')
                    #

                    be, te, step = file_reader(files,
                                               folderpath,
                                               report_name,
                                               column_be,
                                               column_te)

                else:

                    #
                    print(
                        ' -   Picking bindingEnergy as metric')
                    #

                    be, te, step = file_reader(files,
                                               folderpath,
                                               report_name,
                                               column_be,
                                               column_te)

            else:

                be, te, step = file_reader(files,
                                           folderpath,
                                           report_name,
                                           column,
                                           column_te)

        return be, te, step, column_be

    def boltzmann_weighted(be,
                           te,
                           T):
        """
        Function
        ----------
        Calculates boltzmann weighted energy.

        Parameters
        ----------
        - be : list
            Binding energies of all the simulation.
        - te : list
            Total energies of all the simulation.
        - T : float
            Temperature to perform the Boltzmann weights with.
        - steps : list
            Steps associated to poses for all the simulation.

        Returns
        ----------
        - ene_bz : float
            Value of the boltzmann weighted energy.
        """

        exp_bz = np.exp(-te/(R*T))
        nominator = be.dot(exp_bz)
        denominator = np.sum(exp_bz)
        ene_bz = nominator/denominator

        return ene_bz

    def pelesteps_retriever():
        """
        Function
        ----------
        Reads the pelesteps from adaptive.conf or pele.conf.

        Returns
        ----------
        - pele_steps : int
            Number of pele steps of the simulation(s).
        """

        path = str(pathlib.Path().absolute())
        adaptive_path = os.path.join(path, 'adaptive.conf')

        if os.path.isfile(adaptive_path) == True:

            #
            print('     -   Retrieving information from adaptive.conf.')
            #

            with open(adaptive_path) as filein:

                for line in filein:

                    if "peleSteps" in line:

                        peleSteps_string = line.split()[2]
                        pele_steps = int(peleSteps_string.split(',')[0])

        else:

            peleconf_path = os.path.join(path, 'pele.conf')

            #
            print('     -   Retrieving information from pele.conf.')
            #

            if os.path.isfile(peleconf_path) == True:

                with open(peleconf_path) as filein:

                    for line in filein:

                        if "numberOfPeleSteps" in line:

                            pele_steps = int(line.split()[-1])

            else:

                pele_steps = None

                #
                print('     -   No .conf was found.')
                print('     -   The step weighted scoring function will not \n'
                      '         be taken into account.')
                #

        return pele_steps

    def simulation_evolution(report_name,
                             pele_steps,
                             column_be,
                             folderpath):
        """
        Function
        ----------
        Calculates the evolution of the boltzmann weighted value  and the minimum value
        of the metric chosen with column and stores data in two dictionaries.

        Parameters
        ----------
        - files : list
            The path to the directory where the output of the PELE simulation is located.
        - report_name : str
            Name of the reports to obtain the data from
        - pele_steps : int
            Number of pele steps in each epoch.
        - column : int
            Column where the interesting data is located.
        - folderpath : str
            Path where the different epochs of the simulation are located.

        Returns
        ----------
        - bz_dict : dict
            Dictionary with the acummulated boltzmann weighted values of the chosen metric
            per epoch.
        - min_dict : dict
            Dictionary with the accumulated minium values of the chosen metric per epoch.
        - step_dict :dict
            Dictionary with the number of maximum accepted steps per epoch.
        """

        def column_retriever(file):
            """
            Function
            ----------
            Retrieves the position of the binding energy and current energy.

            Parameters
            ----------
            - file : list
                The path to the one report.

            Returns
            ----------
            - column_te : int
                Column where the total energy data is located.
            """

            cont = 0

            with open(file, 'r') as filein:

                for line in filein:

                    if cont == 0:

                        line = line.split('    ')
                        column_te = line.index('currentEnergy') + 1

                    cont += 1

            return column_te

        def step_data_reader(files,
                             report_name,
                             column,
                             column_te,
                             folderpath,
                             step,
                             step_data):
            """
            Function
            ----------
            Opens all reports and extracts information from the column inputed by the user
            at a certain step.

            Parameters
            ----------
            - files : list
                The path to the directory where the output of the PELE simulation is located.
            - report_name : str
                Name of the reports to obtain the data from
            - column : int
                Column where the interesting data is located.
            - column_te : int
                Column where the total energy is located in the report.
            - folderpath : str
                Path where the different epochs of the simulation are located.
            - step : int
                Accepted step we are retrieving data from.
            - step_data : list
                List with all the accumulated data of previous steps.

            Returns
            ----------
            - step_data : list
                List with all the accumulated data of previous steps.
            - step : int
                Accepted step next iteration is going to retrieve data from.
            - te : list
                List with all the total energies obtained.
            """

            for file in files:

                cont = 0

                if file.startswith(report_name):

                    with open(os.path.join(folderpath, file)) as filein:

                        for line in filein:

                            if cont != 0:

                                line = line.split()

                                if step == int(line[2]):

                                    step_data.append(float(line[column - 1]))
                                    te.append(float(line[column_te - 1]))

                            cont += 1

                        else:

                            continue

            step += 1

            return step_data, step, te

        files = os.listdir(folderpath)
        numeric_files = [s for s in os.listdir(folderpath) if s.isnumeric()]

        step_dict = {}
        bz_dict = {}
        min_dict = {}

        if len(numeric_files) != 0:

            column_file = os.path.join(
                folderpath, numeric_files[0], report_name + '_1')
            column_te = column_retriever(column_file)

            #
            print(' -   Calculating...')
            #

            for document in numeric_files:

                step = 0

                step_data = []
                step_list = []
                te = []
                minimum_energy_list = []
                bz_energy_list = []

                new_directory = os.path.join(folderpath, document)

                if os.path.isdir(new_directory) and document.isnumeric():

                    files = os.listdir(new_directory)

                    while step < pele_steps:

                        if step > 1 and bz_energy_list[-1] == bz_energy_list[-2]:

                            break

                        step_list.append(int(step))
                        step_data, step, te = step_data_reader(files,
                                                               report_name,
                                                               column,
                                                               column_te,
                                                               new_directory,
                                                               step,
                                                               step_data)

                        minimum_energy = min(te)

                        ene_bz = boltzmann_weighted(
                            np.array(step_data), np.array(te) - minimum_energy, T)
                        minimum_energy = min(np.array(step_data))

                        bz_energy_list.append(ene_bz)
                        minimum_energy_list.append(minimum_energy)

                step_dict[int(document)] = step_list
                bz_dict[int(document)] = bz_energy_list
                min_dict[int(document)] = minimum_energy_list

        else:

            step = 0

            minimum_energy_list = []
            bz_energy_list = []
            step_data = []
            step_list = []
            te = []

            column_file = os.path.join(
                folderpath, report_name + '_1')
            column_te = column_retriever(column_file)

            #
            print(' -   Calculating...')
            #

            while step < pele_steps:

                if step > 1 and bz_energy_list[-1] == bz_energy_list[-2]:

                    print(
                        ' -   Maximum number of accepted steps in any report:', step - 1)

                    break

                step_list.append(int(step))
                step_data, step, te = step_data_reader(files,
                                                       report_name,
                                                       column_te,
                                                       column_te,
                                                       folderpath,
                                                       step,
                                                       step_data)

                ene_bz = boltzmann_weighted(np.array(step_data), np.array(te), T)
                minimum_energy = min(np.array(step_data))

                bz_energy_list.append(ene_bz)
                minimum_energy_list.append(minimum_energy)

            step_dict[0] = step_list
            bz_dict[0] = bz_energy_list
            min_dict[0] = minimum_energy_list

        return bz_dict, min_dict, step_dict

    def plotter(step_dict,
                bz_dict,
                min_dict,
                path,
                input_folder):
        """
        Function
        ----------
        Orders dictionaries by epochs and plots the data calculated to save the results
        in a newu directory /evolution.

        Parameters
        ----------
        - step_dict :dict
            Dictionary with the number of maximum accepted steps per epoch.
        - bz_dict : dict
            Dictionary with the acummulated boltzmann weighted values of the chosen metric
            per epoch.
        - min_dict : dict
            Dictionary with the accumulated minium values of the chosen metric per epoch.
        - path : str
            Path to the working directory.
        - input_folder : str
            Name of the output directory from which we have obtained the data
            of the dictionaries.
        """

        path_plots = os.path.join(path, 'evolution')

        if os.path.exists(path_plots) == False:
            os.mkdir(path_plots)

        step_dict = dict(sorted(step_dict.items()))
        bz_dict = dict(sorted(bz_dict.items()))
        min_dict = dict(sorted(min_dict.items()))

        plt.figure(figsize=(7, 3.5))

        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()

        for key in step_dict:

            ax1.plot(step_dict[key], bz_dict[key], label=key)
            ax1.set_xlabel('Number of PELEsteps')
            ax1.set_xlabel('Boltzmann weighted metric')
            ax1.set_title('Simultion evolution')
            ax1.legend(loc='best', fontsize=5)
            ax2.plot(step_dict[key], min_dict[key], label=key)
            ax2.legend(loc='best', fontsize=5)

        if len(step_dict[key]) > 100:

            left, bottom, width, height = [0.65, 0.3, 0.2, 0.2]
            ax1_1 = fig1.add_axes([left, bottom, width, height])
            ax1_1.plot(step_dict[key], bz_dict[key], label=key)
            ax1_1.set_xlim([0, 20])

        fig1.savefig(os.path.join(
            path_plots, input_folder + '_bz_evolution.png'))
        fig2.savefig(os.path.join(
            path_plots, input_folder + '_min_evolution.png'))

    def ensambler(input_folder,
                  report_name,
                  column,
                  T,
                  action):
        """
        Function
        ----------
        Function that joins all the other functions.
        """

        folderpath, path = path_definer(input_folder)
        be, te, _, column_be = reader(folderpath,
                                      report_name)

        if action == 'all':

            if column_be != 1:

                min_energy = min(te)
                te = np.array(te) - min_energy

                minimum_energy = min(be)
                be = np.array(be)

                average = np.average(be)

                ene_bz = boltzmann_weighted(be, te, T)

            else:

                minimum_energy = min(te)
                te = np.array(te)

                average = np.average(te)

                ene_bz = boltzmann_weighted(te, te, T)

            #
            print(' ')
            print(' RESULTS:')
            print(' -   Minimum Binding Energy:', minimum_energy)
            print(' -   Average Binding Energy:', average)
            print(' -   Boltzmann weighted Energy:', ene_bz)
            print(' ')
            #
            output_path = os.path.dirname(input_folder)
            output_path = '/'.join(output_path.split('/')[0:-1])
            with open('%s/energy.csv'%output_path, 'w') as fileout:
                fileout.writelines(
                    'Minimum,Average,Boltzmann weighted,Step weighted,Step-Boltzmann weighted,Boltzmann weighted corrected\n'
                    '' + str(minimum_energy) + ',' +
                    str(average) + ',' + str(ene_bz) + '\n'
                )

        elif action == 'evolution':

            #
            print(' -   Retrieving number of pele steps.')
            #

            pele_steps = pelesteps_retriever()

            #
            print(' -   Pele steps:', pele_steps)
            print(' ')
            #

            if column is None:

                bz_dict, min_dict, step_dict = simulation_evolution(report_name,
                                                                    pele_steps,
                                                                    column_be,
                                                                    folderpath)

            else:

                bz_dict, min_dict, step_dict = simulation_evolution(report_name,
                                                                    pele_steps,
                                                                    column,
                                                                    folderpath)

            plotter(step_dict, bz_dict, min_dict, path, input_folder)

            #
            print(' -   Images have been stored in /evolution directory.\n')

        else:

            raise Exception(
                'ActionError: The action you chose is not either all of evolution.')

    #
    print(' ')
    print('**************************************************************')
    print('*                           DiSC                             *')
    print('**************************************************************')
    print(' ')
    #

    ensambler(input_folder,
              report_name,
              column,
              T,
              action)


def main(args):
    """
    Function
    ----------
    It runs statistics function.

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    statistics(input_folder=args.input_folder,
               report_name=args.report_name,
               column=args.column,
               T=args.temperature,
               action=args.action)

    print('*******************************************************************')


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
