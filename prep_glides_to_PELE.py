import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Given a directory\
            with multiple compounds (LIGSdir) docked to a given target\
            (target_pdb) using Glide, prepare them for a PELE simulation.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--LIGSdir',
                                    dest="LIGSdir",
                                    help = "Directory with the docked\
                                            compounds (This must be in PDB format)",
                                    required=True)
    requiredArguments.add_argument('--target_pdb',
                                    dest="target_pdb",
                                    help = "Target PDB", 
                                    required=True)
    requiredArguments.add_argument('--pref',
                                    dest="prefix",
                                    help = "Prefix to save the COMPLEXES which will\
                                            be used as input for the PELE simulation",
                                    required=True)
    requiredArguments.add_argument('--outdir',
                                    dest = 'outdir',
                                    help = 'output directory',
                                    required = True)
    requiredArguments.add_argument('--partition',
                                    dest='partition',
                                    help='MN5 partition either gpp or acc',
                                    required=True)
    parser.add_argument('--HBlist',
                        dest='HBlist',
                        help='Generate a list of all Hydrogen Bonds of the complex',
                        action='store_true',
                        default=False)
    args = parser.parse_args()


    softHB = 'biotite'

    #Parse inputs
    LIGSdir = args.LIGSdir
    target_pdb = args.target_pdb
    prefix = args.prefix
    outdir = args.outdir
    if outdir[-1] != '/': outdir = outdir + '/'
    HBlist = args.HBlist
    partition = args.partition
    if partition == 'gpp':
        qos = 'gp'
    elif partition == 'acc':
        qos = 'acc'
    else:
        raise ValueError('Partition must be either gpp or acc')

    cmd = 'python scripts/prep_files.py --LIGSdir %s --target_pdb %s\
            --pref %s --outdir %s' % (LIGSdir, target_pdb, prefix, outdir)
    if HBlist:
        cmd += ' --HBlist'

    with open('%sprepligs.sh' % outdir, 'w') as fileout:
        fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH --job-name=prepligs\n'
                '#SBATCH --output=%sprepligs.out\n'
                '#SBATCH --error=%sprepligs.err\n'
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
                '%s' % (outdir, outdir, qos, cmd))
    
    os.system('sbatch -A bsc72 %sprepligs.sh' % outdir)
    os.system('rm %sprepligs.*' % outdir)
