import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Given a directory with multiple compounds (LIGSdir) docked to a '
                                                  'given target (target_pdb) using Glide, prepare them for a PELE '
                                                  'simulation.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--LIGSdir', dest="LIGSdir", help = "Directory with the docked compounds (This must be in PDB format)",required=True)
    requiredArguments.add_argument('--target_pdb', dest="target_pdb", help = "Target PDB",required=True)
    requiredArguments.add_argument('-o', dest="outname", help = "Prefix to save the COMPLEXES which will be used as input "
                                                     "for the PELE simulation"
                        ,required=True)
    parser.add_argument('--HBlist',dest='HBlist',help='Generate a list of all Hydrogen Bonds of the complex',action='store_true',default=False)
    requiredArguments.add_argument('--partition', dest='partition', help='MN5 partition either gpp or acc', required=True)

    args = parser.parse_args()


    softHB = 'biotite'
    current_path = os.getcwd()

    #Parse inputs
    LIGSdir = current_path + '/' + args.LIGSdir
    target_pdb = current_path + '/' + args.target_pdb
    complexname = args.outname
    HBlist = args.HBlist
    partition = args.partition
    if partition == 'gpp':
        qos = 'gp'
    elif partition == 'acc':
        qos = 'acc'
    else:
        raise ValueError('Partition must be either gpp or acc')

    cmd = 'python scripts/prep_files.py --LIGSdir %s --target_pdb %s -o %s' % (LIGSdir, target_pdb, complexname)
    if HBlist:
        cmd += ' --HBlist'

    with open('prepligs.sh', 'w') as fileout:
        fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=prepligs\n'
                '#SBATCH --output=prepligs.out\n'
                '#SBATCH --error=prepligs.err\n'
                '#SBATCH --ntasks=1\n'
                '#SBATCH --qos=%s_debug'
                '#SBATCH --time=01:00:00\n'
                '\n'
                'module load anaconda\n'
                'module load intel mkl impi gcc cmake\n'
                'module load transfer\n'
                'module load bsc\n'
                '\n'
                'eval \"$(conda shell.bash hook)\"\n'
                'source activate /gpfs/projects/bsc72/conda_envs/PELE_ConStraParl\n'
                '%s' % (qos, cmd))
    
    os.system('sbatch -A bsc72 prepligs.sh')
    os.system('rm prepligs.*')
