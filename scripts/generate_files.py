import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    parser.add_argument('--LIGSdir', dest="LIGSdir", help = "Directory with the docked compounds",required=True)
    parser.add_argument('-o', dest="outname", help = "",required=True)
    parser.add_argument('--compound', dest="compound", help = "",required=True)
    parser.add_argument('-n', dest="n", help = "Number of processors",required=True)
    parser.add_argument('--HBconsts', dest="HBconsts", help = "HB constrain", nargs='+', action='append', default=None)
    parser.add_argument('--center', dest="center", help = "", nargs='+', default=[None,None,None])
    parser.add_argument('--truncated',dest='truncated',help='',action='store_true',default=False)
    #parser.add_argument('--strain',dest='strain',help='', action='store_true',default=False)
    parser.add_argument('--HBanalysis',dest='HBanalysis',help='', action='store_true',default=False)
    parser.add_argument('--simulation', dest='simulation',help='Choose \'rescoring\' or \'expanded\' simulation type',required=True) 
    args = parser.parse_args()

    #Parse inputs
    ligsdir = args.LIGSdir
    complexesdir = '/'.join(os.path.dirname(ligsdir).split('/')[0:-1]) + '/COMPLEXES'
    ligsdir = '/'.join(os.path.dirname(ligsdir).split('/'))
    outname = args.outname
    compound = args.compound
    n = args.n
    HBconsts = args.HBconsts
    center = args.center
    if len(center) != 3:
        raise ValueError('If given, the center must be 3D')
    truncated = args.truncated
    HBanalysis = args.HBanalysis
    if HBanalysis and not HBconsts:
        raise ValueError('If HBanalysis set as true then you need at least one HBconst')
    simulation = args.simulation

    current_dir = os.getcwd()

    if truncated:
        runinp0 = open('%s/templates/run_template_0_trunc'%current_dir,'r')
    else:
        runinp0 = open('%s/templates/run_template_0'%current_dir,'r')

    runout0 = open('%s/runs/%s/run_%s_0'%(current_dir,outname,compound),'w')
    runinp1 = open('%s/templates/run_template_1'%current_dir,'r')
    runout1 = open('%s/runs/%s/run_%s_1'%(current_dir, outname,compound),'w')

    #Create specific run_0 using run_template_0 as template
    for line in runinp0:
        if '$LIGSDIR' in line:
            line = line.replace('$LIGSDIR',ligsdir)
        if '$OUTNAME' in line:
            line = line.replace('$OUTNAME',outname)
        if '$COMPOUND' in line:
            line = line.replace('$COMPOUND',compound)
        if '$PROCESSORS' in line:
            line = line.replace('$PROCESSORS',n)
        if '$CURRENT':
            line = line.replace('$CURRENT',current_dir)
        runout0.write(line)

    runinp0.close()
    runout0.close()

    #Create specific run_1 using run_template_1 as template
    for line in runinp1:
        if '$LIGSDIR' in line:
            line = line.replace('$LIGSDIR',ligsdir)
        if '$OUTNAME' in line:
            line = line.replace('$OUTNAME',outname)
        if '$COMPOUND' in line:
            line = line.replace('$COMPOUND',compound)
        if '$PROCESSORS' in line:
            line = line.replace('$PROCESSORS',n)
        if '$CURRENT':                                                                                                                                                                             line = line.replace('$CURRENT',current_dir)
        runout1.write(line)

    runinp1.close()
    runout1.close()


    if HBanalysis:
        runout1 = open('%s/runs/%s/run_%s_1'%(current_dir,outname,compound),'a')
        for const in HBconsts:
            runout1.write('\n')
            chain,residue,atom = const[0].split('-')
            cmd = 'python /gpfs/projects/bsc72/COVID/COVID_VS_analysis/FilteringAndClustering.py %s/results/%s/%s/ -n %s --ie_col 5 --rmsd_col 7 -t output/topologies/conntopology_0.pdb -b 2.5 -g2 %s:%s:%s --minimum_g2_conditions 1 -o filtering_results_HB --generate_plots --hbonds_path hbonds.out'%(current_dir,outname,compound,n,chain,residue,atom)
            runout1.write(cmd)
            runout1.close()

    os.system('chmod +x %s/runs/%s/run_%s_0'%(current_dir, outname,compound))
    os.system('chmod +x %s/runs/%s/run_%s_1'%(current_dir, outname,compound))

    #if truncated:
    #    yamlinp = open('%s/templates/yaml_template_trunc.yaml'%current_dir,'r')
   
    if simulation == 'rescoring':
        yamlinp = open('%s/templates/yaml_template.yaml'%current_dir,'r')
    elif simulation == 'expanded':
        yamlinp = open('%s/templates/yaml_template_expanded.yaml'%current_dir,'r')
    else:
        raise ValueError('Indicate the simulation type (rescoring or expanded)')
    yamlout = open('%s/results/%s/%s.yaml'%(current_dir,outname,compound),'w')

    #Create a yaml file
    for line in yamlinp:
        if '$COMPLEXESDIR' in line:
            line = line.replace('$COMPLEXESDIR',complexesdir)
        if '$OUTNAME' in line:
            line = line.replace('$OUTNAME',outname)
        if '$COMPOUND' in line:
            line = line.replace('$COMPOUND',compound)
        if '$PROCESSORS':
            line = line.replace('$PROCESSORS',n)
        if '$CURRENT':
            line = line.replace('$CURRENT',current_dir)
        yamlout.write(line)

    if HBconsts:
        for i,const in enumerate(HBconsts):
            chain0,residue0,atom0 = const[0].split('-')
            chain1,residue1,atom1 = const[1].split('-')
            if i == 0: yamlout.write('atom_dist:\n')
            yamlout.write('  - \"%s:%s:%s\"\n'%(chain0,residue0[3:],atom0))
            yamlout.write('  - \"%s:%s:%s\"\n'%(chain1,residue1[3:],atom1))

    if center!=[None,None,None]:
        yamlout.write('box_center:\n')
        yamlout.write('  - %s\n'%center[0])
        yamlout.write('  - %s\n'%center[1])
        yamlout.write('  - %s\n'%center[2])

    yamlinp.close()
    yamlout.close()
