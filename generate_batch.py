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

    current_dir = os.getcwd()

    #Create runs and results directory
    if not os.path.isdir('%s/runs'%current_dir):
        os.mkdir('%s/runs'%current_dir)
    if not os.path.isdir('%s/runs/%s'%(current_dir, outname)):
        os.mkdir('%s/runs/%s'%(current_dir, outname))
    if not os.path.isdir('%s/results'%current_dir):
        os.mkdir('%s/results'%current_dir)
    if not os.path.isdir('%s/results/%s'%(current_dir, outname)):
        os.mkdir('%s/results/%s'%(current_dir, outname))


    batchfile0 = open('%s/runs/%s_batch_0.sh'%(current_dir, outname),'w')
    batchfile1 = open('%s/runs/%s_batch_1.sh'%(current_dir, outname),'w')
    if strain:
        batchfile2 = open('%s/runs/%s_batch_2.sh'%(current_dir, outname),'w')
        batchfile3 = open('%s/runs/%s_batch_3.sh'%(current_dir, outname),'w')
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
        batchfile0.write('sbatch %s/runs/%s/run_%s_0\n'%(current_dir,outname,compound))
        batchfile1.write('sbatch %s/runs/%s/run_%s_1\n'%(current_dir,outname,compound))
        if strain:
            batchfile2.write('python %s/scripts/ligand_minimization.py -f %s/%s/%s_prep.pdb -d %s/results/%s/%s -r LIG -lf %s/results/%s/%s_min\n'%(current_dir, current_dir, os.path.dirname(ligsdir),compound, current_dir, outname,compound,current_dir, outname,compound))
            batchfile3.write('python %s/scripts/disc.py -d %s/results/%s/%s_min/output/ -c 4\n'%(current_dir,current_dir,outname,compound))
            batchfile3.write('python %s/scripts/corrector.py -d %s/results/%s/%s -lf %s/results/%s/%s_min --skip_strain_per_cluster\n'%(current_dir, current_dir, outname,compound, current_dir, outname,compound))
    batchfile0.close()
    batchfile1.close()
    if strain:
        batchfile2.close()
        batchfile3.close()
