# -*- coding: utf-8 -*-
"""
This module is designed to calculate different scoring methods
for the binding energy of a simulation. Different Scorings Calculator
(DiSC).
"""

# Imports
import sys
import os
import numpy as np
import time
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
from datetime import datetime
import pathlib

# ----------------------------------------------------------------------- #
# Constants:
T = 298.
R = 1.985e-3
K = 1.9872e-3
# ----------------------------------------------------------------------- #


def boltzmann_weighted(df,metric_label,total_energy_label,T=298.):
    """
    Function
    ----------
    Calculates boltzmann weighted ovservable.

    Parameters
    ----------
    - df : pd.DataFrame
        Dataframe with all the information of one simulation protein-ligand.
    - metric_label : str
        Label of the column where the observable you want to measure is located
        (i.e. 'Binding Energy').
    - total_energy_label : str
        Label of the column where the total energy is stored.
    - T : float
        Temperature at which the Boltzmann wheights are going to be calculated.
        Default = 298 K.

    Returns
    ----------
    - ene_bz : float
        Value of the boltzmann weighted energy.
    """

    be = df[metric_label].values
    te = df[total_energy_label].values

    try:
        residency = df['Residence'].values
    except KeyError:
        residency = None
        print('WARNING: Residence not added to the Dataframe. Highly recommendable to use residence_calculator function to have an accurate calculation.')

    # Adding residency
    if residency is not None:
        for be_value, te_value, residency_value in zip(be,te,residency):
            be = np.append(be, np.array((residency_value - 1)*[be_value]))
            te = np.append(te, np.array((residency_value - 1)*[te_value]))

    # Setting 0 of energy at the minimum total energy
    te = te - np.min(te)

    te = te - np.min(te)
    exp_bz = np.exp(-te/(R*T))
    nominator = be.dot(exp_bz)
    denominator = np.sum(exp_bz)
    ene_bz = nominator/denominator

    return ene_bz

def bindingFreeEnergy(df,metric_label='BindingEnergy',total_energy_label='currentEnergy',T=298.):
    total_energy = df[total_energy_label]
    relative_energy = total_energy-total_energy.min()
    Z = np.sum(np.exp(-relative_energy/(K*T)))
    probability = np.exp(-relative_energy/(K*T))/Z
    bfe = np.sum(probability*df[metric_label])

    return bfe

def clustering(weight,
               angles,
               rotatable_bonds,
               path_results,
               path_images):
    """
    Function
    ----------
    Cluster the results obtained and stored in a data frame.

    Parameters
    ----------
    - n_cluster : int
        Number of clusters to cluster data.
    - clustering_method : str
        Method to cluster data: bin or kmeans.
    - simulation_df : pd.DataFrame
        Data frame with all the rotatable bonds' dihedral angle values of all the simulation
        with corresponding model, trajectory and epoch.
    - dihedral_bond_df : pd.DataFrame
        Data frame with rotatable bonds, atoms conforming it and the index assigned.
    - path_results : str
        Path to the directory where the results will be stored.

    """

    def binning(weight,
                angles,
                rotatable_bonds,
                path_results,
                path_images):
        """
        Function
        ----------
        Cluster the results obtained from simulation by binning the results.
        Entropic contributions are calculated from the binned data.

        Parameters
        ----------
        - simulation_df : pd.DataFrame
            Data frame with all the rotatable bonds' dihedral angle values of all the simulation
            with corresponding model, trajectory and epoch.
        - path_images : str
            Path to the directory where the results will be stored.
        """

        entropy_contribution = []

        results = defaultdict(list)

        for rot_bond, value, weight in zip(rotatable_bonds, angles, weight):
            results[rot_bond].append([value,weight])

        rot_bond_values = list(results.items())

        for rot_bond, values in rot_bond_values:

            print('     -   Dihedral angle being clustered {dihedral} out of {num_of_dihedrals}'.format(dihedral=rot_bond, num_of_dihedrals=max(rotatable_bonds)), end="\r")

            angles_vector = np.array(np.array(values)[:,0])
            weights_vector = np.array(np.array(values)[:,1])

            bin_edges = np.histogram_bin_edges(
                angles_vector, bins=10)
            density, _ = np.histogram(
                angles_vector, bins=bin_edges, density=True, weights=weights_vector)
            dense_bins = density[density != 0]

            entropy_contribution.append(
                np.sum(np.array([p*np.log(p) for p in dense_bins])))

            # Plot
            plt.title('Dihedral ' + str(rot_bond) + ' distribution')
            plt.hist(values,
                     bins=bin_edges, density=True)
            plt.xlabel('Dihedral angle (º)')
            plt.ylabel('Density')
            plt.xlim(-180, 180)
            plt.xticks(list(np.arange(-180, 190, 30)))
            plt.savefig(os.path.join(path_images, 'dihedral_' +
                                     str(rot_bond) + '_distribution.png'), format='png', transparent=True)
            plt.close()

        entropy_contributions = np.array(entropy_contribution)
        S = -(R)*np.sum(entropy_contributions)

        entropy_percentages = 100. * \
            np.array(entropy_contribution)/np.sum(entropy_contributions)
        entropy_df = pd.DataFrame(entropy_percentages)
        entropy_df = entropy_df.round(decimals=2).T
        entropy_df.columns = [
            'Dihedral_' + str(value) + '_%' for value in list(np.arange(1, len(rot_bond_values) + 1))]
        entropy_df.insert(loc=0, column='S (kcal/mol K)',
                          value="{:3.6f}".format(S))
        entropy_df.to_csv(os.path.join(
            path_results, 'entropy.csv'), index=False)

        print(' -   Entropic information written in /dihedrals/entropy.csv.')

    binning(weight,
            angles,
            rotatable_bonds,
            path_results,
            path_images)

def bootstrapping(number_of_samples,
                  metric,
                  original_df,
                  snapshots):
    """
    Function
    ----------
    Bootstrap function that performs all the calculations.

    Parameters
    ----------
    - input_folder : str
        Name of the output directory where the simulation
        is located.
    - report_name : str
        Name of the report files used for the simulation.
    - number_of_samples : int
        Number of bootstrap datasets to generate.
    - metric : str
        Name of the metric you are interested in.
    """

    def path_definer():
        """
        Function
        ----------
        Defines the important paths that are going to be used in this function.

        Parameters
        ----------

        Returns
        ----------
        - path_results : str
            Path to the results folder of the analysis.
        """

        path = str(pathlib.Path().absolute())
        path_results = os.path.join(path, 'analyzer','bootstrap')

        if os.path.isdir(path_results) is False:
            os.mkdir(path_results)

        return path_results

    def bootstrap_dataset(original_df,
                          snapshots):
        """
        Function
        ----------
        Generates a bootstrap datset from the data retrieved and random numbers.

        Parameters
        ----------
        - original_df : pd.DataFrame
            Data frame with all the data from the reports.
        - snapshots : int
            Number of total snapshots in a simulation.

        Return
        ----------
        - bootstrap_df : pd.DataFrame
            Data frame with all the data of the bootstrap dataset.
        """

        def random_numbers_generator(snapshots):
            """
            Function
            ----------
            Generates a random array of length: snapshots, of numbers
            ranging from 0 to: snapshots.

            Parameters
            ----------
            - snapshots : int
                Number of snaphots in the simulation.

            Return
            ----------
            - random_numbers : np.array
                Array of length: snapshots, of numbers ranging [0,snapshots].
            """

            seed = int(''.join(c for c in str(
                datetime.now().time()) if c.isdigit())[-5:])
            np.random.seed(seed)
            random_numbers = np.random.randint(0, snapshots, snapshots)

            return random_numbers

        def random_to_data(original_df,
                           random_numbers):
            """
            Function
            ----------
            Generates a dataframe from the original dataframe and the
            random values of the array.

            Parameters
            ----------
            - original_df : pd.DataFrame
                Data frame with all the data from the reports.
            - random_numbers : np.array
                Array of length: snapshots, of numbers ranging [0,snapshots].

            Return
            ----------
            - bootstrap_df : pd.DataFrame
                Data frame with all the data from original dataframe and the
                random numbers.
            """

            indexes_list = [[number*6, number*6 + 1, number*6 + 2, number *
                             6 + 3, number*6 + 4, number*6 + 5] for number in random_numbers]
            indexes = [item for sublist in sorted(
                indexes_list) for item in sublist]
            bootstrap_df = original_df.iloc[indexes]

            return bootstrap_df

        random_numbers = random_numbers_generator(snapshots)
        bootstrap_df = random_to_data(original_df,
                                      random_numbers)

        return bootstrap_df

    def metric_calculator(metric,
                          bootstrap_df):
        """
        Function
        ----------
        Calculates all the scores for a given dataframe and metric.

        Parameters
        ----------
        - metric : str
            Metric the user is interested in.
        - bootstrap_df : pd.DataFrame
            Data frame with all the data from original dataframe and the
            random numbers.

        Return
        ----------
        - minimum : float
            The dataframe's minimum value for the specific metric.
        - maximum : float
            The dataframe's maximum value for the specific metric.
        - average : float
            The dataframe's average value for the specific metric.
        - boltzmann : float
            The dataframe's boltzmann average value for the specific metric.
        """

        def from_dataframe_to_vector(metric,
                                     bootstrap_df):
            """
            Function
            ----------
            From the dataframe obtains the values of the interesting metric and
            the total energy.

            Parameters
            ----------
            - metric : str
                Metric the user is interested in.
            - bootstrap_df : pd.DataFrame
                Data frame with all the data from original dataframe and the
                random numbers.

            Return
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the
                metric of interest.
            - total_ene_vector : np.array
                Array with all the values of the dataframe corresponding to the
                total energy.
            """

            metric_vector = bootstrap_df[bootstrap_df['metric']
                                         == metric]['value'].to_numpy()
            total_ene_vector = bootstrap_df[bootstrap_df['metric']
                                            == 'currentEnergy']['value'].to_numpy()

            return metric_vector, total_ene_vector

        def minimum_score(metric_vector):
            """
            Function
            ----------
            Compute minimum of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the
                metric of interest.

            Return
            ----------
            - minimum : float
                Minimum value of the inputed vector.
            """

            minimum = min(list(metric_vector))

            return minimum

        def maximum_score(metric_vector):
            """
            Function
            ----------
            Compute maximum of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the
                metric of interest.

            Return
            ----------
            - maximum : float
                Maximum value of the inputed vector.
            """

            maximum = max(list(metric_vector))

            return maximum

        def average_score(metric_vector):
            """
            Function
            ----------
            Compute average of a vector.

            Parameters
            ----------
            - metric_vector: np.array
                Array with all the values of the dataframe corresponding to the
                metric of interest.

            Return
            ----------
            - average : float
                Average value of the inputed vector.
            """

            average = np.average(metric_vector)

            return average

        def boltzmann_score(metric_vector, total_ene_vector, T):
            """
            Function
            ----------
            Compute boltzmann average of a vector.

            Parameters
            ----------
            - metric_vector : np.array
                Array with all the values of the dataframe corresponding to the
                metric of interest.
            - T : float
                Temperature to perform the Boltzmann weights with.
            - total_ene_vector : np.array
                Steps associated to poses for all the simulation.

            Returns
            ----------
            - boltzmann : float
                Boltzmann average value of the inputed vectors.
            """

            total_ene_vector -= min(list(total_ene_vector))

            exp_bz = np.exp(-total_ene_vector/(R*T))
            nominator = metric_vector.dot(exp_bz)
            denominator = np.sum(exp_bz)
            boltzmann = nominator/denominator

            return boltzmann

        metric_vector, total_ene_vector = from_dataframe_to_vector(
            metric, bootstrap_df)

        minimum = minimum_score(metric_vector)
        maximum = maximum_score(metric_vector)
        average = average_score(metric_vector)
        boltzmann = boltzmann_score(metric_vector, total_ene_vector, T)

        return minimum, maximum, average, boltzmann

    def iterator(number_of_samples,
                 metric,
                 original_df,
                 snapshots):
        """
        Function
        ----------
        Iterates number_of_samples times to obtain the different scores for
        the specified metric.

        Parameters
        ----------
        - number_of_samples : int
            Number of times we want to iterate the process of creating
            a bootsrap dataset and calculate scores.
        - metric : str
            Metric the user is interested in.
        - original_df : pd.DataFrame
            Data frame with all the data from the reports.
        - snapshots : int
            Number of total snapshots in a simulation.

        Return
        ----------
        - average : float
            The dataframe's minimum value for the specific metric.
        - maximum : float
            The dataframe's maximum value for the specific metric.
        - average : float
            The dataframe's average value for the specific metric.
        - boltzmann : float
            The dataframe's boltzmann average value for the specific metric.
        """

        minimum_list = []
        maximum_list = []
        average_list = []
        boltzmann_list = []

        for i in range(number_of_samples):

            bootstrap_df = bootstrap_dataset(original_df, snapshots)
            minimum, maximum, average, boltzmann = metric_calculator(
                metric, bootstrap_df)

            minimum_list.append(minimum)
            maximum_list.append(maximum)
            average_list.append(average)
            boltzmann_list.append(boltzmann)

        return minimum_list, maximum_list, average_list, boltzmann_list

    def statistics(metric_list):
        """
        Function
        ----------
        Calculates average and error of the inputed list.

        Parameters
        ----------
        - metric_list : list
            List of values we want to know the average and the standard
            deviation of.

        Return
        ----------
        - average : float
            Avergage of the inputed list.
        - standard_dev : float
            Standard deviation of the inputed list.
        """

        vector = np.array(metric_list)
        average = np.average(vector)
        standard_dev = np.std(vector)

        return average, standard_dev

    def plotter(path_results,
                data,
                average,
                error,
                scoring):
        """
        Function
        ----------
        Plots results obtained with bootstrap.

        Parameters
        ----------
        - path_results : str
            Path to the results folder of the bootstrap analysis.
        - data : list
            List of data points calculated from bootstrap we want to plot.
        - average : float
            Avergage of the inputed list.
        - error : float
            Standard deviation of the inputed list.
        - scoring : string
            Name of the scoring function used.
        """

        from seaborn import displot

        displot(data, kind="kde", color='black', label='KDE plot')
        plt.title(scoring + ' Distribution')
        plt.axvline(x=average, color='red', label='Average = ' +
                    str("{:.3f}".format(average)))
        plt.axvline(x=average + error, color='green',
                    label='Error = ' + str("{:.3f}".format(error)))
        plt.axvline(x=average - error, color='green')
        plt.legend(loc="best")
        plt.xlabel(scoring)
        plt.ylabel('Density')
        plt.tight_layout()
        plt.savefig(os.path.join(path_results, scoring +
                    '_distribution.png'), format='png', bbox_inches="tight")

    path_results = path_definer()

    #
    print('\n *   Bootstrapping')
    print(' -   Generating ' + str(number_of_samples) + ' datasets...')
    #

    minimum_list,\
        maximum_list,\
        average_list,\
        boltzmann_list = iterator(number_of_samples,
                                  metric,
                                  original_df,
                                  snapshots)

    average_minimum, standard_dev_minimum = statistics(minimum_list)
    average_maximum, standard_dev_maximum = statistics(maximum_list)
    average_average, standard_dev_average = statistics(average_list)
    average_boltzmann, standard_dev_boltzmann = statistics(boltzmann_list)

    data = {
        'minimum': [average_minimum, standard_dev_minimum],
        'maximum': [average_maximum, standard_dev_maximum],
        'average': [average_average, standard_dev_average],
        'boltzmann': [average_boltzmann, standard_dev_boltzmann],
    }

    #
    print(' -   Writing files')
    #

    results_df = pd.DataFrame(data, index=['average', 'error'])
    results_df.to_csv(os.path.join(path_results, 'results.csv'))

    plotter(path_results,
            minimum_list,
            average_minimum,
            standard_dev_minimum,
            'Minimum')

    plotter(path_results,
            maximum_list,
            average_maximum,
            standard_dev_maximum,
            'Maximum')

    plotter(path_results,
            average_list,
            average_average,
            standard_dev_average,
            'Average')

    plotter(path_results,
            boltzmann_list,
            average_boltzmann,
            standard_dev_boltzmann,
            'Boltzmann')

def residence_calculator(df,pele_steps):
    """
    Function
    ----------
    Computes residency of the models and adds a column with residency.

    Parameters
    ----------
    - df : pd.DataFrame
        Dataframe with all the information of one simulation protein-ligand.
    - pele_steps : int
        Number of pele stepes per epoch of the simulation.

    Returns
    ----------
    - df_residence : pd.DataFrame
        Dataframe with all the information of one simulation protein-ligand.
    """

    # Retrieveing accepted steps
    accepted_pele_steps = np.array(df.Step.to_numpy())

    auxiliar = np.array(accepted_pele_steps[1:])
    auxiliar[auxiliar == 0] = pele_steps
    following_model_accepted_pele_steps = np.append(auxiliar,pele_steps)

    # Calculating residence
    residence = following_model_accepted_pele_steps - accepted_pele_steps

    # Adding residence
    df_residence = df.copy()
    df_residence['Residence'] = residence

    return df_residence






























































