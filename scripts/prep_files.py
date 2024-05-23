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
    #Parse inputs
    LIGSdir = args.LIGSdir
    if LIGSdir[-1] != '/': LIGSdir = LIGSdir + '/'
    target_pdb = args.target_pdb
    complexname = args.outname
    data = args.csv
    feature_filter = args.feature_filter
    if args.value_filter != None:
        value_filter = float(args.value_filter)
    sdf_out = args.sdf_out
    HBlist = args.HBlist

    TARGdir = os.path.dirname(target_pdb)
    structureTARG = parsePDB('%s'%(target_pdb))

    #If missing create a COMPLEX directory
    COMPLEXESdir = '/'.join(os.path.dirname(LIGSdir).split('/')[0:-1]) + '/COMPLEXES'
    if not os.path.isdir(COMPLEXESdir):
        os.mkdir(COMPLEXESdir)
    LIGS = glob.glob('%s/*.pdb'%LIGSdir)

    #If needed create an HBlist directory
    if HBlist:
        HBdir = '/'.join(os.path.dirname(LIGSdir).split('/')[0:-1]) + '/HBlists'
    if not os.path.isdir(HBdir):
        os.mkdir(HBdir)

    #Apply Glide feature filter by creating an auxiliary dataframe (datafiltered)
    if feature_filter != None:
        data = sts.loadSchrodingerCSV(data)
        data = data[['Job Name','Title',feature_filter]]
        datafiltered = data[(data[feature_filter] <= value_filter )]
        print(datafiltered)
        print('- Filter of %.3f applied on %s '%(value_filter,feature_filter))
        print('- From %d without filter to %d compounds with filter'%(len(set(data['Title'].tolist())),len(set(datafiltered['Title'].tolist()))))
    else:
        datafiltered = data

    count_preparedLIGS = 0
    prepLIGs = []
    for LIG in LIGS:
        name = os.path.basename(LIG)
        name = name.replace('.pdb','')
        name = name.split('_') #To remove the identifier of multiple docking poses/chirality for the same LIG
        fullname = '_'.join(name)
        if len(name) > 1:
            name = '_'.join(name[0:-1])
        else:
            name = '_'.join(name)
        if feature_filter != None:
            dataaux = datafiltered[(datafiltered['Title'] == name)] #Get the rows representing all compound poses
            if dataaux.shape[0] == 0:
                if feature_filter != None:
                    print('%s didn\'t pass the %s filter'%(name,feature_filter))
                else:
                    print('%s wasn\'t in the maestro data table'%name)
                continue
            else:
                idx = dataaux[feature_filter].idxmin() #Select the min feature value for a lig with multiple poses/chirality
                print('%s has a %s of %.3f (below %.3f). We prepare it for a further PELE simulation.'%(fullname,feature_filter,datafiltered[feature_filter][idx],value_filter))

        count_preparedLIGS+=1
        #Change LIGS chain name to L and resid to LIG
        structureLIG = parsePDB('%s'%(LIG))
        hvLIG = structureLIG.getHierView()
        for chain in hvLIG:
            chain.setChid('L')
            for res in chain:
                res.setResname('LIG')
        prepLIG = '%s%s_prep.pdb'%(LIGSdir,fullname)
        writePDB(prepLIG,structureLIG)
        prepLIGs.append(prepLIG)

        #Merge Target and Ligand
        structureCOMPLEX = structureTARG + structureLIG
        prepCOMPLEX = '%s/%s_%s.pdb'%(COMPLEXESdir,complexname,fullname)
        writePDB(prepCOMPLEX,structureCOMPLEX)

        #If asked compute a list of all Hydrogen Bonds of the complex
        if HBlist:
            if softHB == 'biotite':
                pdb_file = pdb.PDBFile.read(prepCOMPLEX)
                pdbcomplex = pdb_file.get_structure()
                triplets, mask = hbonds(pdbcomplex)
                triplets = np.squeeze(triplets)
                f = open('%s/%s.txt'%(HBdir,fullname),'w')
                for hbond in triplets:
                    D_idx = hbond[0]
                    A_idx = hbond[2]
                    D = pdbcomplex[0,D_idx].res_name + str(pdbcomplex[0,D_idx].res_id) + '-' + pdbcomplex[0,D_idx].atom_name
                    A = pdbcomplex[0,A_idx].res_name + str(pdbcomplex[0,A_idx].res_id) + '-' +  pdbcomplex[0,A_idx].atom_name
                    f.write('%s -- %s\n'%(D,A))
                f.close()
            elif softHB == 'mdtraj':
                pdbcomplex = md.load_pdb(prepCOMPLEX)
                #The criterion employed to compute an hbond is: theta > 120 and d(HAcceptor) < 2.5A both in at least 10% of the trajectory
                hbonds = md.baker_hubbard(pdbcomplex, periodic=False)
                #hbonds = np.squeeze(md.wernet_nilsson(pdbcomplex, periodic=False))
                #List of hbonds (row of three elements: index of the donor atom, index of the hydrogen atom and index of the acceptor atom)
                #Define label function to parse the hbonds
                label = lambda hbond : '%s -- %s' % (pdbcomplex.topology.atom(hbond[0]), pdbcomplex.topology.atom(hbond[2]))
                f = open('%s/%s.txt'%(HBdir,fullname),'w')
                for hbond in hbonds:
                    l = label(hbond)
                    f.write(l+'\n')
                f.close()


        #Add connections to LIGs_prep
        toextend = open('%s/%s_prep.pdb'%(LIGSdir,fullname),'a')
        toinclude = open('%s/%s.pdb'%(LIGSdir,fullname),'r')
        for line in toinclude:
            if 'CONECT' in line or 'END' in line:
                toextend.write(line)
        toextend.close()
        toinclude.close()

    if sdf_out != None:
        DB = MolDB(pdbList=prepLIGs,chirality=False)
        print('%d prepared molecules loaded into a MolDB object'%len(DB.dicDB.keys()))
        DB.filter_similarity()
        print('%d prepared molecules after filtering them by similarity'%len(DB.dicDB.keys()))
        for k in DB.dicDB.keys():
            print(DB.dicDB[k][1])
        DB.save_MolDB_sdf(sdf_out)

    print('%s molecules preapred for PELE run'%count_preparedLIGS)
