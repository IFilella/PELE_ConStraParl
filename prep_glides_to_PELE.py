from prody import *
import glob
import os
import argparse
from MolecularAnalysis.analysis import sts
from MolecularAnalysis.moldb import MolDB
import mdtraj as md
import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure.hbond as hbonds

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Given a directory with multiple compounds (LIGSdir) docked to a '
                                                  'given target (target_pdb) using Glide, prepare them for a PELE '
                                                  'simulation. Docked compounds can be filtered by a Glide feature'
                                                  ' (feature_filter and value_filter) when a Glide table is provided'
                                                  ' (csv).')
    parser.add_argument('--LIGSdir', dest="LIGSdir", help = "Directory with the docked compounds (This must be in PDB format)",required=True)
    parser.add_argument('--target_pdb', dest="target_pdb", help = "Target PDB",required=True)
    parser.add_argument('-o', dest="outname", help = "Prefix to save the COMPLEXES which will be used as input "
                                                     "for the PELE simulation"
                        ,required=True)
    parser.add_argument('--HBlist',dest='HBlist',help='Generate a list of all Hydrogen Bonds of the complex',action='store_true',default=False)
    parser.add_argument('--csv',dest="csv",help="Glide feature table",default=None)
    parser.add_argument('--feature_filter',dest="feature_filter",help="Glide feature to apply a compound filter"
                        ,default=None)
    parser.add_argument('--value_filter',dest="value_filter",help="Value to which apply the glide feature filter"
                        ,default=None)
    parser.add_argument('--sdf_out',dest="sdf_out",help="Output sdf file with the molecules fullfiling the Glide feature"
                                                        " filter"
                       ,default=None)

    args = parser.parse_args()


    softHB = 'biotite'
    current_path = os.getcwd()

    #Parse inputs
    LIGSdir = current_path + '/' + args.LIGSdir
    target_pdb = current_path + '/' + args.target_pdb
    complexname = args.outname
    if args.csv:
        data = current_path + '/' + args.csv
    else:
        data = args.csv
    feature_filter = args.feature_filter
    if args.value_filter:
        value_filter = float(args.value_filter)
    else:
        value_filter = args.value_filter
    if args.sdf_out:
        sdf_out = current_path + '/' + args.sdf_out
    else:
        sdf_out = args.sdf_out
    HBlist = args.HBlist

    cmd = 'python scripts/prep_files.py --LIGSdir %s --target_pdb %s -o %s' % (LIGSdir, target_pdb, complexname)
    if HBlist:
        cmd += ' --HBlist'
    if data:
        cmd += ' --csv %s' % data
    if value_filter:
        cmd += ' --value_filter %s' % value_filter
    if feature_filter:
        cmd += ' --feature_filter %s' % feature_filter
    if sdf_out:
        cmd += ' --sdf_out %s' % sdf_out

    with open('prepligs.sh', 'w') as fileout:
        fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=prepligs\n'
                '#SBATCH --output=prepligs.out\n'
                '#SBATCH --error=prepligs.err\n'
                '#SBATCH --ntasks=1\n'
                '#SBATCH --qos=gp_debug'
                '#SBATCH --time=01:00:00\n'
                '\n'
                'module load anaconda\n'
                'module load intel mkl impi gcc cmake\n'
                'module load transfer\n'
                'module load bsc\n'
                '\n'
                'eval \"$(conda shell.bash hook)\"\n'
                'source activate PELE_ConStraParl\n'
                '%s' % cmd)
    
    os.system('sbatch -A bsc72 prepligs.sh')
    os.system('rm prepligs.*')
