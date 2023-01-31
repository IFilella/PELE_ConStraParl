import glob
import argparse
import os
import subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    parser.add_argument('--LIGSdir', dest="LIGSdir", help = "Directory with the docked compounds",required=True)
    parser.add_argument('-o', dest="outname", help = "",required=True)
    parser.add_argument('-n', dest="n", help = "Number of processors",required=True)
    parser.add_argument('--center', dest="center", help = "PELE box center", nargs='+', default=[None,None,None])
    parser.add_argument('--HBfilter', dest="HBfilter", help = "HB filter in the format C-RES000-AA."
                        "Looks for an HydrogenBond with the specified residue",default=None)
    parser.add_argument('--truncated',dest='truncated',help='',action='store_true',default=False)
    parser.add_argument('--strain',dest='strain',help='Create an extra batch_2 file to get the minimal energy of each ligand with PELE', action='store_true',default=False)
    args = parser.parse_args()

    #Parse inputs
    ligsdir = args.LIGSdir
    outname = args.outname
    n = args.n
    HBfilter = args.HBfilter
    if HBfilter != None:
        HBlistdir = '/'.join(os.path.dirname(ligsdir).split('/')[0:-1]) + '/HBlists'
        if not os.path.isdir(HBlistdir):
            raise ValueError('Can\'t pass an HBfilter without an HBlists directory')
    center = args.center
    if len(center) != 3:
        raise ValueError('If given, the center must be 3D')
    truncated = args.truncated
    strain = args.strain

    compounds = glob.glob('%s/*_prep.pdb'%(ligsdir))
    compounds = [os.path.basename(compound).split('_prep.pdb')[0] for compound in compounds]

    #Create runs and results directory
    if not os.path.isdir('runs'):
        os.mkdir('runs')
    if not os.path.isdir('runs/%s'%outname):
        os.mkdir('runs/%s'%outname)
    if not os.path.isdir('results'):
        os.mkdir('results')
    if not os.path.isdir('results/%s'%(outname)):
        os.mkdir('results/%s'%(outname))


    batchfile0 = open('runs/%s_batch_0.sh'%(outname),'w')
    batchfile1 = open('runs/%s_batch_1.sh'%(outname),'w')
    if strain:
        batchfile2 = open('runs/%s_batch_2.sh'%(outname),'w')
        batchfile3 = open('runs/%s_batch_3.sh'%(outname),'w')
    for compound in compounds:
        cmd = 'python scripts/generate_files.py --LIGSdir %s -o %s --compound %s -n %s'%(ligsdir,outname,compound,n)
        if center != [None,None,None]:
            cmd += ' --center %s %s %s'%(center[0],center[1],center[2])
        if HBfilter != None:
            _HBfilter = '-'.join(HBfilter.split('-')[1:])
            chain = HBfilter.split('-')[0]
            grep_process = subprocess.Popen(['grep',_HBfilter,HBlistdir+'/%s.txt'%compound], stdout= subprocess.PIPE)
            grep_out, grep_err = grep_process.communicate()
            grep_out = grep_out.decode("utf-8").replace("\n","")
            if grep_out == "":
                print('Compound %s is not fulfilling the HBfilter %s'%(compound,HBfilter))
                continue
            grep_out = grep_out.split(' -- ')
            _grep_out = ['','']
            for i,element in enumerate(grep_out):
                if 'LIG' in element:  _grep_out[1] = 'L-' + element
                else: _grep_out[0] = chain + '-' + element
            print(_grep_out)
            if _grep_out[1] == '':
                print('Compound %s is not fulfilling the HBfilter %s with LIG'%(compound,HBfilter))
                continue
            cmd += ' --HBfilter %s %s'%(_grep_out[0], _grep_out[1])
        if truncated:
            cmd += ' --truncated'
        #if strain:
        #    cmd += ' --strain'
        print(cmd)
        os.system(cmd)
        batchfile0.write('sbatch %s/run_%s_0\n'%(outname,compound))
        batchfile1.write('sbatch %s/run_%s_1\n'%(outname,compound))
        if strain:
            batchfile2.write('python ../scripts/ligand_minimization.py -f ../%s/%s_prep.pdb -d ../results/%s/%s -r LIG -lf ../results/%s/%s_min\n'%(os.path.dirname(ligsdir),compound,outname,compound,outname,compound))
            batchfile3.write('python ../scripts/disc.py -d ../results/%s/%s_min/output/ -c 4\n'%(outname,compound))
            batchfile3.write('python ../scripts/corrector.py -d ../results/%s/%s -lf ../results/%s/%s_min --skip_strain_per_cluster\n'%(outname,compound,outname,compound))
    batchfile0.close()
    batchfile1.close()
    if strain:
        batchfile2.close()
        batchfile3.close()
