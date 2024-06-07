import glob
import argparse
import os
import subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='This script creates\
            the execution files needed for each ligand to complete\
            the PELE simulation and puts them into a new runs directory.')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--LIGSdir',
                                    dest="LIGSdir",
                                    help = "Directory with the docked compounds",
                                    required=True)
    requiredArguments.add_argument('--pref',
                                    dest="prefix",
                                    help = "should be the same as the one used in\
                                            prep_glide_to_PELE.py", 
                                    required=True)
    requiredArguments.add_argument('--outdir',
                                    dest="outdir",
                                    help = "should be the same as the one used in\
                                            prep_glide_to_PELE.py", 
                                    required=True)
    requiredArguments.add_argument('-n',
                                    dest="n",
                                    help = "Number of processors", 
                                    required=True)
    # parser.add_argument('--center',
    #                     dest="center",
    #                     help = "PELE box center",
    #                     nargs='+',
    #                     default=[None,None,None])
    parser.add_argument('--HBconsts',
                        dest="HBconsts",
                        help = "HB consts in the format C-RES000-AA C-RES000-AA ...\
                                Looks for HydrogenBonds with the specified atoms",
                        default=None,
                        nargs='+')
    parser.add_argument('--strain', 
                        dest='strain', 
                        help='Create an extra batch_2 file to get the minimal energy\
                                of each ligand with PELE',
                        action='store_true', 
                        default=False)
    requiredArguments.add_argument('--partition',
                                    dest='partition',
                                    help='MN5 partition either gpp or acc',
                                    required=True)
    args = parser.parse_args()

    # Parse inputs
    repodir = os.getcwd()
    ligsdir = args.LIGSdir
    os.chdir(ligsdir)
    ligsdir = os.getcwd()
    os.chdir(repodir)
    prefix = args.prefix
    outdir = args.outdir
    os.chdir(outdir)
    outdir = os.getcwd()
    os.chdir(repodir)
    n = args.n
    HBconsts = args.HBconsts
    if HBconsts:
        HBlistdir = outdir + '/HBlists'
        if not os.path.isdir(HBlistdir):
            raise ValueError('Can\'t pass an HBconsts without an HBlists directory\
                    on the outdir directory')
    # center = args.center
    # if len(center) != 3:
    #     raise ValueError('If given, the center must be 3D')
    strain = args.strain
    partition = args.partition
    if partition != 'gpp' and partition!= 'acc':
        raise ValueError('Partition must be either gpp or acc')

    compounds = glob.glob('%s/*_prep.pdb'%(ligsdir))
    compounds = [os.path.basename(compound).split('_prep.pdb')[0] for compound in compounds]

    # Create runs and results directory
    if not os.path.isdir('%s/runs' % outdir):
        os.mkdir('%s/runs' % outdir)
    if not os.path.isdir('%s/results' % outdir):
        os.mkdir('%s/results' % outdir)

    batchfile1 = open('%s/%s_batch_1.sh' % (outdir, prefix), 'w')
    if strain:
        batchfile2 = open('%s/%s_batch_2.sh'%(outdir, prefix), 'w')
        batchfile3 = open('%s/%s_batch_3.sh'%(outdir, prefix), 'w')
    for i, compound in enumerate(compounds):
        print('%d %s'%(i+1, compound))
        cmd = 'python scripts/generate_files.py --LIGSdir %s --pref %s --outdir %s --compound %s -n %s --partition %s' % (ligsdir, prefix, outdir, compound, n, partition)

        # if center != [None,None,None]:
        #     cmd += ' --center %s %s %s' % (center[0], center[1], center[2])

        if HBconsts:
            for const in HBconsts:
                _const = '-'.join(const.split('-')[1:])
                chain = const.split('-')[0]
                grep_process = subprocess.Popen(['grep',
                    _const, HBlistdir + '/%s.txt' %compound], stdout= subprocess.PIPE)
                grep_out, grep_err = grep_process.communicate()
                if grep_out == "":
                    print('Compound %s is not fulfilling the HBconst %s' % (compound, const))
                    continue
                grep_out = grep_out.decode("utf-8").split("\n")
                for aux in grep_out:
                    if aux == "": continue
                    _grep_out = ['', '']
                    aux = aux.split(' -- ')
                    for i, element in enumerate(aux):
                        if 'LIG' in element: _grep_out[1] = 'L-' + element
                        else: _grep_out[0] = chain + '-' + element
                    if _grep_out[1] == '':
                        print('Compound %s is not fulfilling the HBconst %s with\
                                LIG ' % (compound, const))
                    if _grep_out[0] != "" and _grep_out[1] != "":
                        print('Compound %s is going to be const at atom %s to \
                                atom %s' % (compound, _grep_out[1], _grep_out[0]))
                        cmd += ' --HBconsts %s %s' % (_grep_out[0], _grep_out[1])
        if strain:
            cmd += ' --strain'
        os.system(cmd)
        batchfile1.write('sbatch -A bsc72 %s/runs/run_%s_1\n' % (outdir, compound))
        if strain:
            batchfile2.write('sbatch -A bsc72 %s/runs/run_%s_2\n' % (outdir, compound))
            batchfile3.write('sbatch -A bsc72 %s/runs/run_%s_3\n' % (outdir, compound))
        print('-----------------------')
    batchfile1.close()
    if strain:
        batchfile2.close()
        batchfile3.close()
