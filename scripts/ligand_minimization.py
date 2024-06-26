# -*- coding: utf-8 -*-
"""
This module is designed to run PELE wth a single ligand to then cluster
positions and obtain energies.
"""

__author__ = "Ignasi Puch-Giner"
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import os
import pathlib
import argparse
import shutil
from distutils.dir_util import copy_tree


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

    parser.add_argument("-f", "--input_file", type=str, dest="input_file",
                        default=None, required=True, help="Original pdb name of lonely ligand.")
    parser.add_argument("-d", "--directory", type=str, dest="input_folder",
                        default='LIG_Pele', help="Name of the directory where the simulation\
                        is located.")
    parser.add_argument("-r", "--residue_name", type=str, dest="residue_name",
                        default='LIG', help="Ligand's residue name.")
    parser.add_argument("-ff", "--force_field", type=str, dest="force_field",
                        default=None, help="Force field with which the simulation will run.")
    parser.add_argument("-sm", "--solvent_model", type=str, dest="solvent_model",
                        default=None, help="Solvent model with which the simulation will run.")
    parser.add_argument("-lf", "--ligand_folder", type=str, dest="ligand_folder",
                        default='LIG_ligand', help="Name of the output folder of the ligand simulation.")

    parser.add_argument("--skip_simulation", dest="simulation_bool",
                        default=False, action='store_true', help="Flag to choose if the PELE simulation is performed automatically or not.")

    parsed_args = parser.parse_args(args)

    return parsed_args


def ligand_prepare(input_folder,
                  residue_name,
                  pdb_name,
                  force_field,
                  solvent_model,
                  simulation_bool,
                  ligand_folder):
    """
    Function
    ----------
    Prepares the PELE simulation of the ligand.

    Parameters
    ----------
    - input_folder : str
        The path to the directory created by the induced fit simulation.
    - residue_name : str
        Residue name of the ligand in the pdb of each cluster.
    - pdb_name : str
        Name of the pdb we ant to perform the simulation with.
    - force_field : str
        Force field with which the simulation will run.
    - solvent_model : str
        Solvent model with which the simulation will run.
    """

    def path_definer(input_folder,
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
            Residue name of the ligand in the pdb of each cluster.

        Returns
        ----------
        - path : str
            Absolute path from where the code is being executed.
        - path_previous_simulation: str
            The path to the directory generated by the simulation we want to analyze the clusters
            from.
        - path_ligand : str
            The path to the generated directory containing all the necessary files to perform the
            PELE energy calculation.
        """

        path = str(pathlib.Path().absolute())
        path_previous_simulation = os.path.join(path, input_folder)

        if os.path.isdir(path_previous_simulation) == False:
            raise Exception('PathError: There is no folder with this name: '
                            + path_previous_simulation + '. Please check the path and the folder name.')

        path_ligand = os.path.join(path, ligand_folder)

        if os.path.exists(path_ligand) == False:
            os.mkdir(path_ligand)

        return path, path_previous_simulation, path_ligand

    def ff_sm_checker(force_field,
                      solvent_model):
        """
        Function
        ----------
        Checks the value introduced with the flags and if none is given, it gives one.

        Parameters
        ----------
        - force_field : str
            Force field value given by the user.
        - solvent_model : str
            Solvent model value given by the user.

        Returns
        ----------
        - force_field : str
            Force field with which the simulation will run.
        - solvent_model : str
            Solvent model with which the simulation will run.
        """

        if force_field is None:
            print(
                '                              WARNING:                               \n'
                '   No force field was introduced with the flag -ff (i.e. -ff         \n'
                '   OPLS2005).The script automatically sets the forcefield to         \n'
                '   OPLS2005.'
                '\n'
            )
            force_field = 'OPLS2005'

        if solvent_model is None:
            print(
                '                              WARNING:                               \n'
                '   No solvent model was introduced with the flag -sm (i.e. -sm       \n'
                '   OBC). The script automatically sets the solvent model to          \n'
                '   VDGBNP.'
                '\n'
            )
            solvent_model = 'VDGBNP'

        return force_field, solvent_model

    def directory_preparation(path,
                              pdb_name,
                              residue_name,
                              path_ligand,
                              path_previous_simulation):
        """
        Function
        ----------
        Prepares the PELE simulation of the ligand.

        Parameters
        ----------
        - path : str
            Absolute path from where the code is being executed.
        - pdb_name : str
            Name of the pdb we want to perform the simulation with.
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster.
        - path_ligand : str
            The path to the generated directory containing all the necessary files to perform the
            PELE energy calculation.
        - path_previous_simulation: str
            The path to the directory generated by the simulation we want to analyze the clusters
            from.
        """

        new_path = os.path.join(path_ligand, 'ligand.pdb')
        shutil.copy(os.path.join(path, pdb_name), new_path)

        path_DataLocal = os.path.join(path_ligand, 'DataLocal')

        if os.path.exists(path_DataLocal) == False:
            os.mkdir(path_DataLocal)

        path_previous_DataLocal = os.path.join(path_previous_simulation,
                                               'DataLocal')

        if os.path.exists(path_previous_DataLocal) == False:

            print('\n'
                  '                              WARNING:                               \n'
                  '   No DataLocal directory has been found at ' + input_folder + '.    \n'
                  '   The DataLocal folder should be copied at ../' +
                  residue_name + '_ligand/' + '\n'
                  )
        else:

            copy_tree(path_previous_DataLocal, path_DataLocal)

    def constraint_retriever(pdb_name,
                             path,
                             residue_name):
        """
        Function
        ----------
        Retrieves information from pdb to introduce in the .conf file.

        Parameters
        ----------
        - pdb_name : str
            Name of the pdb we want to perform the simulation with.
        - path : str
            Absolute path from where the code is being executed.
        - residue_name : str
            Residue name of the ligand in the pdb of each cluster.


        Returns
        ----------
        - chain : str
            Chain name of the ligand.
        - number : str
            Number of the ligand.
        - atom : str
            First atom of te pdb.
        """

        cont = 0

        with open(os.path.join(path, pdb_name)) as filein:
            for line in filein:
                if cont != 0:
                    if residue_name in line:

                        line = line.split()

                        if len(line[4]) == 5:

                            chain = list(line[4])[0]
                            number = "".join(list(line[4])[1:])
                            #atom = '_' + line[2] + '_'
                            atom = line[2]
                            if len(atom) == 1 or len(atom) == 2:
                                atom = '_' + atom + '_'
                            elif len(atom) == 3:
                                atom = '_' + atom

                        else:

                            chain = line[4]
                            number = line[5]
                            #atom = '_' + line[2] + '_'
                            atom = line[2]
                            if len(atom) == 1 or len(atom) == 2:
                                atom = '_' + atom + '_'
                            elif len(atom) == 3:
                                atom = '_' + atom

                        #if len(atom) == 3:

                        #    atom = atom + '_'

                        break

                cont += 1

        return chain, number, atom

    def write_files(force_field,
                    solvent_model,
                    chain,
                    number,
                    atom,
                    path_ligand,
                    pdb_name):
        """
        Function
        ----------
        Writes all the necessary files to run a PELE simulation and a
        platform analysis.

        Parameters
        ----------
        - force_field : str
            Force field to be used in the upcoming simulation.
        - solvent_model : str
            Solvent model to be used in the upcoming simulation.
        - path_ligand : str
            Path to the folder where the files will be stored.
        """

        path_to_run = os.path.join(path_ligand, 'run')

        with open(os.path.join(path_ligand, 'pele.conf'), 'w') as fileout:

            fileout.writelines(
                '{\n'
                '  "licenseDirectoryPath": "/gpfs/projects/bsc72/PELE++/license",\n'
                '  "Initialization": {\n'
                '    "Complex": {\n'
                '      "files": [\n'
                '        {\n'
                '          "path": "' + path_ligand + '/ligand.pdb"\n'
                '        }\n'
                '      ]\n'
                '    },\n'
                '    "ForceField": "' + force_field + '",\n'
                '    "Solvent": {\n'
                '      "ionicStrength": 0.15,\n'
                '      "solventType": "' + solvent_model + '",\n'
                '      "useDebyeLength": true\n'
                '    }\n'
                '  },\n'
                '  "verboseMode": false,\n'
                '  "commands": [\n'
                '    {\n'
                '      "commandType": "peleSimulation",\n'
                '      "RandomGenerator": {\n'
                '        "seed": 12345\n'
                '      },\n'
                '      "selectionToPerturb": {\n'
                '        "chains": {\n'
                '          "names": [\n'
                '            "L"\n'
                '          ]\n'
                '        }\n'
                '      },\n'
                '      "PELE_Output": {\n'
                '        "savingFrequencyForAcceptedSteps": 1,\n'
                '        "savingMode": "savingTrajectory",\n'
                '        "reportPath": "output/report",\n'
                '        "trajectoryPath": "output/trajectory.pdb"\n'
                '      },\n'
                '      "PELE_Parameters": {\n'
                '        "anmFrequency": 0,\n'
                '        "sideChainPredictionFrequency": 2,\n'
                '        "minimizationFrequency": 1,\n'
                '        "waterPerturbationFrequency": 1,\n'
                '        "perturbationCOMConstraintConstant": 0,\n'
                '        "sideChainPredictionRegionRadius": 6,\n'
                '        "activateProximityDetection": true,\n'
                '        "temperature": 1500,\n'
                '        "numberOfPeleSteps": 2000\n'
                '      },\n'
                '\n'
                '      "constraints":[\n'
                '      { "type": "constrainAtomToPosition", "springConstant": 0.5, "equilibriumDistance": 0.0, "constrainThisAtom": "' +
                chain + ':' + number + ':' + atom + '" }\n'
                '      ],\n'
                '\n'
                '      "SideChainPerturbation": {\n'
                '        "sideChainsToPerturb": {\n'
                '          "links": {\n'
                '            "ids": [\n'
                '              "' + chain + ':' + number + '"\n'
                '            ]\n'
                '          }\n'
                '        },\n'
                '        "parameters": {\n'
                '          "overlapFactor": 0.5,\n'
                '          "numberOfTrials": 20,\n'
                '          "atLeastOneSelectedTrial": true,\n'
                '          "maxTrialsForAtLeastOne": 50\n'
                '        }\n'
                '      },\n'
                '      "ANM": {\n'
                '        "algorithm": "NULL"\n'
                '        },\n'
                '\n'
                '      "Minimizer": {\n'
                '        "algorithm": "TruncatedNewton",\n'
                '        "parameters": {\n'
                '          "MinimumRMS": 0.2,\n'
                '          "alphaUpdated": false,\n'
                '          "nonBondingListUpdatedEachMinStep": true\n'
                '        }\n'
                '      },\n'
                '      "PeleTasks": [\n'
                '        {\n'
                '          "metrics": [\n'
                '            {\n'
                '              "type": "sasa",\n'
                '              "tag": "sasaLig",\n'
                '              "selection": {\n'
                '                "chains": {\n'
                '                  "names": [\n'
                '                    "L"\n'
                '                  ]\n'
                '                }\n'
                '              }\n'
                '            }\n'
                '          ]\n'
                '        }\n'
                '      ]\n'
                '    }\n'
                '  ]\n'
                '}\n'
            )

        with open(path_to_run, 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH -J ' + pdb_name.split('.pdb')[0] + '\n'
                '#SBATCH --output=PELEne.out\n'
                '#SBATCH --error=PELEne.err\n'
                '#SBATCH --ntasks=48\n'
                '#SBATCH --qos=gp_bscls\n'
                '#SBATCH --time=02:00:00\n'
                '\n'
                'module purge\n'
                'module load anaconda\n'
                'module load intel mkl impi gcc cmake\n'
                'module load transfer\n'
                'module load bsc\n'
                '\n'
                'mpirun -np 48 /gpfs/projects/bsc72/PELE++/mnv/1.8.1b1/bin/PELE_mpi --control-file ' +
                path_ligand + '/pele.conf --license-directory /gpfs/projects/bsc72/PELE++/license\n'
            )

        with open(os.path.join(path_ligand, 'run_analysis'), 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=analysis\n'
                '#SBATCH --output=analysis.out\n'
                '#SBATCH --error=analysis.err\n'
                '#SBATCH --ntasks=48\n'
                '#SBATCH --time=00-05:00:00\n'
                '\n'
                'module load anaconda\n'
                'module load intel mkl impi gcc cmake\n'
                'module load transfer\n'
                'module load bsc\n'
                '\n'
                'eval "$(conda shell.bash hook)"\n'
                'conda activate /gpfs/projects/bsc72/conda_envs/platform\n'
                '\n'
                'python script.py\n'
            )

        with open(os.path.join(path_ligand, 'script.py'), 'w') as fileout:

            fileout.writelines(
                'from pele_platform.analysis import Analysis\n'
                '\n'
                'analysis = Analysis(resname="' + residue_name +
                '", chain="L", simulation_output="output", report="report", cpus=48)\n'
                'analysis.generate(path="analysis", clustering_type="meanshift")\n'
            )

        return path_to_run

    #
    print(' ')
    print('*******************************************************************')
    print('*                           peleLInEn                              *')
    print('* --------------------------------------------------------------- *')
    print('*      Ligand\'s internal energy from ligand PELE simulation       *')
    print('*******************************************************************')
    print(' ')
    #

    path, path_previous_simulation, path_ligand = \
        path_definer(input_folder, residue_name)

    force_field, solvent_model = ff_sm_checker(force_field, solvent_model)

    print(' -   Writing necessary files to run a PELE simulation and')
    print('     a PELE platform analysis afterwards.')
    print(' -   Copying DataLocal and generating ligand.pdb file.\n')

    directory_preparation(path,
                          pdb_name,
                          residue_name,
                          path_ligand,
                          path_previous_simulation)

    chain, number, atom = constraint_retriever(pdb_name, path, residue_name)

    path_to_run = write_files(force_field,
                              solvent_model,
                              chain,
                              number,
                              atom,
                              path_ligand,
                              pdb_name)

    if not simulation_bool:

        print(' ')

        os.chdir(path_ligand)
        os.system("sbatch -A bsc72 %s" % path_to_run)

    else:

        print(' -   PELE simulation skipped.')

    #
    print(' ')


def main(args):
    """
    Function
    ----------
    It reads the command-line arguments and runs ligand_results.
    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    ligand_prepare(input_folder=args.input_folder,
                  residue_name=args.residue_name,
                  pdb_name=args.input_file,
                  force_field=args.force_field,
                  solvent_model=args.solvent_model,
                  simulation_bool=args.simulation_bool,
                  ligand_folder=args.ligand_folder)


if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)
