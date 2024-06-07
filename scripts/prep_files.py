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
    parser = argparse.ArgumentParser(description ='Given a directory\
            with multiple compounds (LIGSdir) docked to a given target\
            (target_pdb) using Glide, prepare them for a PELE simulation.')
    parser.add_argument('--LIGSdir',
                        dest="LIGSdir",
                        help = "Directory with the docked compounds\
                                (This must be in PDB format)",
                        required=True)
    parser.add_argument('--target_pdb',
                        dest="target_pdb",
                        help = "Target PDB",
                        required=True)
    parser.add_argument('--pref',
                        dest="prefix",
                        help = "Prefix to save the COMPLEXES which will\
                                be used as input for the PELE simulation", 
                                required=True)
    parser.add_argument('--outdir',
                        dest = 'outdir',
                        help = 'output directory',
                        required = True)
    parser.add_argument('--HBlist', 
                        dest='HBlist', 
                        help='Generate a list of all Hydrogen Bonds of the complex', 
                        action='store_true', 
                        default=False)
    args = parser.parse_args()

    # Parse inputs
    LIGSdir = args.LIGSdir
    if LIGSdir[-1] != '/': LIGSdir = LIGSdir + '/'
    target_pdb = args.target_pdb
    prefix = args.prefix
    outdir = args.outdir
    if outdir[-1] != '/': outdir = outdir + '/'
    HBlist = args.HBlist

    TARGdir = os.path.dirname(target_pdb)
    structureTARG = parsePDB('%s'%(target_pdb))

    # If missing create a COMPLEX directory
    COMPLEXESdir = outdir + 'COMPLEXES'
    
    if not os.path.isdir(COMPLEXESdir):
        os.mkdir(COMPLEXESdir)
    LIGS = glob.glob('%s/*.pdb'%LIGSdir)
    
    # If needed create an HBlist directory
    if HBlist:
        HBdir = outdir + 'HBlists'
        if not os.path.isdir(HBdir):
            os.mkdir(HBdir)
    
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

        count_preparedLIGS += 1
        # Change LIGS chain name to L and resid to LIG
        structureLIG = parsePDB('%s'%(LIG))
        hvLIG = structureLIG.getHierView()
        for chain in hvLIG:
            chain.setChid('L')
            for res in chain:
                res.setResname('LIG')
        prepLIG = '%s%s_prep.pdb'%(LIGSdir, fullname)
        writePDB(prepLIG, structureLIG)
        prepLIGs.append(prepLIG)

        #Merge Target and Ligand
        structureCOMPLEX = structureTARG + structureLIG
        prepCOMPLEX = '%s/%s_%s.pdb' % (COMPLEXESdir, prefix, fullname)
        writePDB(prepCOMPLEX, structureCOMPLEX)

        #If asked compute a list of all Hydrogen Bonds of the complex
        if HBlist:
            pdb_file = pdb.PDBFile.read(prepCOMPLEX)
            pdbcomplex = pdb_file.get_structure()
            triplets, mask = hbonds(pdbcomplex)
            triplets = np.squeeze(triplets)
            f = open('%s/%s.txt'%(HBdir, fullname),'w')
            for hbond in triplets:
                D_idx = hbond[0]
                A_idx = hbond[2]
                D = pdbcomplex[0,D_idx].res_name + str(pdbcomplex[0,D_idx].res_id) + '-' + pdbcomplex[0,D_idx].atom_name
                A = pdbcomplex[0,A_idx].res_name + str(pdbcomplex[0,A_idx].res_id) + '-' +  pdbcomplex[0,A_idx].atom_name
                f.write('%s -- %s\n'%(D,A))
            f.close()


        #Add connections to LIGs_prep
        toextend = open('%s/%s_prep.pdb'%(LIGSdir, fullname),'a')
        toinclude = open('%s/%s.pdb'%(LIGSdir, fullname),'r')
        for line in toinclude:
            if 'CONECT' in line or 'END' in line:
                toextend.write(line)
        toextend.close()
        toinclude.close()

    print('%s molecules preapred for PELE run'%count_preparedLIGS)
