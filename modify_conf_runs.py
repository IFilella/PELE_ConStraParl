import json
import argparse
import os
import glob
import subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    parser.add_argument('--RESdir', dest="resdir", help = "Directory with the simulations input files prior to the PELE simulation",required=True)
    parser.add_argument('--LIGSdir', dest="LIGSdir", help = "Directory with the docked compounds",required=True)
    parser.add_argument('--strain', dest="strain", help = 'To add internal energy of the ligand in the report needed to compute the strain energy', default=False, action = 'store_true')
    parser.add_argument('--HBfilter', dest="HBfilter", help = "HB filter",default=None)
    parser.add_argument('--sprC', dest="springConstant", help = "springConstant in case of HBfilter adding a constraint",default=5)
    parser.add_argument('--eqD', dest="equilibriumDistance", help = "equilibriumDistance in case of HBfilter adding a constraint",default=2.7)
    args = parser.parse_args()

    #Parse
    resdir = args.resdir
    HBfilter = args.HBfilter
    ligsdir = args.LIGSdir
    if HBfilter != None:
        HBlistdir = '/'.join(os.path.dirname(ligsdir).split('/')[0:-1]) + '/HBlists'
        print(HBlistdir)
        if not os.path.isdir(HBlistdir):
            raise ValueError('Can\'t pass an HBfilter without an HBlists directory')
    springConstant = args.springConstant
    equilibriumDistance = args.equilibriumDistance
    strain = args.strain
    
    current_dir = os.getcwd()    

    compounds = glob.glob('%s/%s/*.yaml'%(current_dir,resdir))
    compounds = [os.path.basename(comp).replace('.yaml','') for comp in compounds]

    #Change debug mode to false in yaml files
    for compound in compounds:
        with open(current_dir + '/' + resdir+compound+".yaml", "r") as sources:
            lines = sources.readlines()
        with open(current_dir + '/' + resdir+compound+".yaml", "w") as sources:
            for line in lines:
                sources.write(line.replace('debug: true','restart: true'))

    #If HBfilter then add a constraint in pele.conf
    if HBfilter != None:
        chain0,res0,atom0 = HBfilter.split('-')
        for compound in compounds:
            grep_process = subprocess.Popen(['grep',res0+'-'+atom0,current_dir + '/' + HBlistdir+'/%s.txt'%compound], stdout=subprocess.PIPE)
            grep_out, grep_err = grep_process.communicate()
            grep_out = grep_out.decode("utf-8").replace("\n","")
            if grep_out == "":
                raise ValueError('Compound %s is not fulfilling the HBfilter %s'%(compound,HBfilter))
            grep_out = grep_out.split(' -- ')
            for i,element in enumerate(grep_out):
                if 'LIG' in element:
                    chain1 = 'L'
                    element = element.split('-')
                    res1 = element[0]
                    atom1 = element[1]
            atomlist = [atom0,atom1]
            for i,atom in enumerate(atomlist):
                if len(atom) == 1:
                    atomlist[i] = '_'+atom+'__'
                elif len(atom) == 2:
                    atomlist[i] = '_'+atom+'_'
                elif len(atom) == 3:
                    atomlist[i] = '_'+atom
                elif len(atom) == 4:
                    atomlist[i] = atom
            _atom0, _atom1 = atomlist
            const = '{ \"type\": \"constrainAtomsDistance\", \"springConstant\": %s, \"equilibriumDistance\": %s, \"constrainThisAtom\": \"%s:%s:%s\", \"toThisOtherAtom\": \"%s:%s:%s\"  },'%(springConstant,equilibriumDistance,chain1,res1[3:],_atom1,chain0,res0[3:],_atom0)
            print(compound)
            peleconf = current_dir + '/' + resdir+compound+'/pele.conf'
            with open(peleconf, "r") as infile:
                lines = infile.readlines()
            with open(peleconf, "w") as outfile:
                for line in lines:
                    if '\"constraints\"' in line:
                        line += '\t\t\t' + const + '\n'
                    outfile.write(line)

    if strain:
        metric = '\n                        { \"type\": \"internalEnergy\",\n                           \"atomSetSelection\": { \"chains\": { \"names\": [\"L\"]  }  }\n                        },\n\n'
        for compound in compounds:
            print(compound)
            peleconf = current_dir + '/' + resdir+compound+'/pele.conf'
            with open(peleconf, "r") as infile:
                lines = infile.readlines()
            with open(peleconf, "w") as outfile:
                for line in lines:
                    if '\"metrics\"' in line:
                        line+=metric
                    outfile.write(line)
