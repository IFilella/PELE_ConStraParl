import argparse
import os

def to_PELE_selection_format(chain, residue, atom):
    if len(atom) == 1:
        _atom = '_'+atom+'__'
    elif len(atom) == 2:
        _atom = '_'+atom+'_'
    elif len(atom) == 3:
        _atom = '_'+atom
    elif len(atom) == 4:
        _atom = atom
    _residue = residue[3:]
    return '%s:%s:%s' % (chain, _residue, _atom)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    requiredArguments = parser.add_argument_group('required arguments')
    requiredArguments.add_argument('--LIGSdir',
                                    dest="LIGSdir", 
                                    help = "Directory with the docked compounds", 
                                    required=True)
    requiredArguments.add_argument('--pref',
                                    dest="prefix",
                                    help = "",
                                    required=True)
    requiredArguments.add_argument('--outdir',
                                    dest="outdir",
                                    help = "",
                                    required=True)
    requiredArguments.add_argument('--compound',
                                    dest="compound",
                                    help = "", 
                                    required=True)
    requiredArguments.add_argument('-n',
                                    dest="n",
                                    help = "Number of processors", 
                                    required=True)
    parser.add_argument('--HBconsts',
                        dest="HBconsts",
                        help = "HB constrain",
                        nargs='+',
                        action='append',
                        default=None)
    # parser.add_argument('--center',
    #                     dest="center",
    #                     help = "",
    #                     nargs='+',
    #                     default=[None,None,None])
    parser.add_argument('--strain',
                        dest='strain',
                        help='',
                        action='store_true', 
                        default=False)
    requiredArguments.add_argument('--partition',
                                    dest='partition',
                                    help='MN5 partition either gpp or acc',
                                    required=True)
    args = parser.parse_args()

    #Parse inputs
    ligsdir = args.LIGSdir
    prefix = args.prefix
    outdir = args.outdir
    complexesdir = outdir + '/COMPLEXES'
    compound = args.compound
    n = args.n
    HBconsts = args.HBconsts
    # center = args.center
    # if len(center) != 3:
    #     raise ValueError('If given, the center must be 3D')
    strain = args.strain
    partition = args.partition
    if partition == 'gpp':
        qos = 'gp'
    elif partition == 'acc':
        qos = 'acc'
    else:
        raise ValueError('Partition must be either gpp or acc')

    repodir = os.getcwd()

    runinp1 = open('%s/templates/run_template_1' % repodir, 'r')
    runout1 = open('%s/runs/run_%s_1' % (outdir, compound), 'w')

    #Create specific run_1 using run_template_1 as template
    for line in runinp1:
        if '$LIGSDIR' in line:
            line = line.replace('$LIGSDIR',ligsdir)
        if '$PREFIX' in line:
            line = line.replace('$PREFIX', prefix)
        if '$COMPOUND' in line:
            line = line.replace('$COMPOUND',compound)
        if '$PROCESSORS' in line:
            line = line.replace('$PROCESSORS',n)
        if '$OUTDIR':
            line = line.replace('$OUTDIR', outdir)
        if '$QOS' in line:
            line = line.replace('$QOS', qos)
        if '$REPO' in line:
            line = line.replace('$REPO', repodir)
        runout1.write(line)

    os.system('chmod +x %s/runs/run_%s_1'%(outdir, compound))
    runinp1.close()
    runout1.close()

    if strain:
        runinp2 = open('%s/templates/run_template_2' % repodir, 'r')
        runout2 = open('%s/runs/run_%s_2'%(outdir, compound),'w')

        for line in runinp2:
            if '$LIGSDIR' in line:
                line = line.replace('$LIGSDIR',ligsdir)
            if '$PREFIX' in line:
                line = line.replace('$PREFIX', prefix)
            if '$COMPOUND' in line:
                line = line.replace('$COMPOUND', compound)
            if '$OUTDIR':
                line = line.replace('$OUTDIR', outdir)
            if '$QOS' in line:
                line = line.replace('$QOS', qos)
            if '$REPO' in line:
                line = line.replace('$REPO', repodir)
            runout2.write(line)
    
        os.system('chmod +x %s/runs/run_%s_2'%(outdir, compound))
        runinp2.close()
        runout2.close()

        runinp3 = open('%s/templates/run_template_3' % repodir,'r')
        runout3 = open('%s/runs/run_%s_3' % (outdir, compound),'w')
 
        for line in runinp3:
            if '$LIGSDIR' in line:
                line = line.replace('$LIGSDIR', ligsdir)
            if '$PREFIX' in line:
                line = line.replace('$PREFIX', prefix)
            if '$OUTDIR':
                line = line.replace('$OUTDIR', outdir)
            if '$COMPOUND' in line:
                line = line.replace('$COMPOUND', compound)
            if '$PROCESSORS' in line:
                line = line.replace('$PROCESSORS', n)
            if '$REPO':
                line = line.replace('$REPO', repodir)
            if '$QOS' in line:
                line = line.replace('$QOS', qos)
            runout3.write(line)
    
        os.system('chmod +x %s/runs/run_%s_3'%(outdir, compound))
        runinp3.close()
        runout3.close()
 
    yamlinp = open('%s/templates/yaml_template.yaml' % repodir,'r')
    yamlout = open('%s/results/%s.yaml' % (outdir, compound),'w')

    #Create a yaml file
    for line in yamlinp:
        if '$COMPLEXESDIR' in line:
            line = line.replace('$COMPLEXESDIR', complexesdir)
        if '$PREFIX' in line:
            line = line.replace('$PREFIX', prefix)
        if '$OUTDIR':
            line = line.replace('$OUTDIR', outdir)
        if '$COMPOUND' in line:
            line = line.replace('$COMPOUND', compound)
        if '$PROCESSORS':
            line = line.replace('$PROCESSORS', n)
        if '$REPO':
            line = line.replace('$REPO', repodir)
        if 'adaptive_epochs' in line:
            if strain or HBconsts:
                line += '    pele_tasks_metrics:\n'
            if strain:
                line += '      - type: internal_energy\n'
                line += '        tag: InternalEnergy\n'
                line += '        atomset_selection: L\n'
            if HBconsts:
                for i, const in enumerate(HBconsts): 
                    chain0, residue0, atom0 = const[0].split('-')
                    chain1, residue1, atom1 = const[1].split('-')
                    sel0 = to_PELE_selection_format(chain0, residue0, atom0)
                    sel1 = to_PELE_selection_format(chain1, residue1, atom1)
                    line += '      - type: com_distance\n'
                    line += '        tag: %s%s%s%s%s%s\n' % (chain1, residue1[3:], atom1, chain0, residue0[3:], atom0)
                    line += '        selection_group_1: \"%s\"\n' % sel1
                    line += '        selection_group_2: \"%s\"\n' % sel0
                line += '    pele_harmonic_constraints:\n'
                for i, const in enumerate(HBconsts):    
                    chain0, residue0, atom0 = const[0].split('-')
                    chain1, residue1, atom1 = const[1].split('-')
                    sel0 = to_PELE_selection_format(chain0, residue0, atom0)
                    sel1 = to_PELE_selection_format(chain1, residue1, atom1)
                    line += '      - type: atom_atom\n'
                    line += '        spring_constant: 5.0\n'
                    line += '        equilibrium_distance: 2.7\n'
                    line += '        atom_selection1: \"%s\"\n' % sel1
                    line += '        atom_selection2: \"%s\"\n' % sel0  
        yamlout.write(line)

    # if center!=[None,None,None]:
    #     yamlout.write('box_center:\n')
    #     yamlout.write('  - %s\n'%center[0])
    #     yamlout.write('  - %s\n'%center[1])
    #     yamlout.write('  - %s\n'%center[2])

    yamlinp.close()
    yamlout.close()
