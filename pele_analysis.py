from MolecularAnalysis.mol import Mol
import argparse
import glob
import pandas as pd
import numpy as np
import os
import sys
sys.path.append('scripts/')
import functions
from rdkit import Chem

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    parser.add_argument('--RESdir', dest="RESdir", help = "All PELE results directory", required=True)
    parser.add_argument('--LIGSdir',dest="LIGSdir",help = 'Original directory with Ligands',required=True)
    parser.add_argument('-o', dest="outname", help = "Outname",required=True)
    parser.add_argument('--glide',dest='glide',help='Glide docking csv file', default=None)
    parser.add_argument('--PELEmetric',dest='metric',help='Final PELE metric',default='all_BFE')
    parser.add_argument('--charge',dest='charge',help='if added compute the net charge of each molecule (only works if glide.csv with SMILES is provided)', default=False,action='store_true')
    parser.add_argument('--maxtanimoto',dest='maxtanimoto',help='add it if maxtanimoto property included in the glide.csv', default=False,action='store_true')
    parser.add_argument('--batch',dest='batch',help='add it if batch property included in the glide.csv', default=False,action='store_true')
    parser.add_argument('--smile',dest='smile',help='add it if SMILE property included in the glide.csv', default=False,action='store_true')

    args = parser.parse_args()

    pele_output = args.RESdir
    if pele_output[-1] != '/': pele_output += '/'
    ligs_dir = args.LIGSdir
    if ligs_dir[-1] != '/': ligs_dir += '/'
    outname = args.outname
    glide_csv = args.glide
    metric = args.metric
    charge = args.charge
    maxtanimoto = args.maxtanimoto
    batch = args.batch
    smile = args.smile

    molecules_out = glob.glob('%s/*/'%pele_output)
    molecules_out = [molecule_out for molecule_out in molecules_out if '_min' not in molecule_out]
    molecules = [molecule.split('/')[-2] for molecule in molecules_out]

    #Get a .csv of all PELE simulations with PELE metrics and if required glide gscore
    f_out =  open('%s_pelemetrics.csv'%outname,'w')
    header = 'molecule,smile,batch,maxtanimoto,charge,glidegscore,all_avgBE,all_minBE,all_BZM,all_BFE,pop_avgBE,pop_minBE,pop_BZM,pop_BFE,lowavg_avgBE,lowavg_minBE,lowavg_BZM,lowavg_BFE,lowBE_avgBE,lowBE_minBE,lowBE_BZM,lowBE_BFE,lowBZM_avgBE,lowBZM_minBE,lowBZM_BZM,lowBZM_BFE,lowBFE_avgBE,lowBFE_minBE,lowBFE_BZM,lowBFE_BFE'
    if glide_csv:
        if not smile:
            header = header.replace('smile,','')
        if not charge:
            header = header.replace('charge,','')
        if not batch:
            header = header.replace('batch,','')
        if not maxtanimoto:
            header = header.replace('maxtanimoto,','')
        glide_df = pd.read_csv(glide_csv)
    else:
        header = header.replace('smile,','').replace('charge,','').replace('batch,','').replace('maxtanimoto,','').replace('glidegscore,','')
    f_out.write(header+'\n')
    for i,molecule in enumerate(molecules):
        print('-', i+1, molecule)
        #Glide gscore
        if glide_csv:
            aux = molecule.split('_')
            if len(aux) == 4:
                row = glide_df[(glide_df['Title'] == '_'.join(aux[0:3]))]
            else:
                row = glide_df[(glide_df['Title'] == molecule)]
            glidegscore = row['glide gscore'].min()
            if maxtanimoto:
                maxt = row['maxtanimoto'].max()
            if smile:
                sml = row['smiles'].iloc[0]
            if batch:
                bat = row['batch'].iloc[0]
            if charge:
                mol = Mol(smile=smile)
                netcharge = Chem.rdmolops.GetFormalCharge(mol.mol)
        #PELE data
        pele_data = '%sanalysis/data.csv'%(molecules_out[i])
        try:
            pele_df = pd.read_csv(pele_data)
        except:
            print('PELE ERROR')
            continue

        #Metrics for the complete PELE simulation
        pele_df = functions.residence_calculator(pele_df,12)
        all_BFE = functions.bindingFreeEnergy(pele_df)
        all_avgBE, all_minBE, all_BZM = functions.calculateMetrics(pele_df)

        #Get metrics for most populated PELE simulation cluster
        pop_df, pop_avgBE, pop_minBE, pop_BZM = functions.cluster_pop(pele_df)
        pop_BFE = functions.bindingFreeEnergy(pop_df)

        lowavg_dict, lowBE_dict, lowBZM_dict, lowBFE_dict = functions.cluster_energy(pele_df)

        lowavg_dict_key = list(lowavg_dict.keys())[0]
        lowavg_avgBE, lowavg_minBE, lowavg_BZM, lowavg_BFE = lowavg_dict[lowavg_dict_key]

        lowBE_dict_key = list(lowBE_dict.keys())[0]
        lowBE_avgBE, lowBE_minBE, lowBE_BZM, lowBE_BFE = lowBE_dict[lowBE_dict_key]

        lowBZM_dict_key = list(lowBZM_dict.keys())[0]
        lowBZM_avgBE, lowBZM_minBE, lowBZM_BZM, lowBZM_BFE = lowBZM_dict[lowBZM_dict_key]

        lowBFE_dict_key = list(lowBFE_dict.keys())[0]
        lowBFE_avgBE, lowBFE_minBE, lowBFE_BZM, lowBFE_BFE = lowBFE_dict[lowBFE_dict_key]

        f_out.write('%s,'%molecule)
        if glide_csv:
            if smile: f_out.write('%s,'%sml)
            if batch: f_out.write('%s,'%bat)
            if maxtanimoto: f_out.write('%.4f,'%maxt)
            if charge: f_out.write('%d,'%netcharge)
            f_out.write('%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n'%(glidegscore,all_avgBE,all_minBE,all_BZM,all_BFE,pop_avgBE,pop_minBE,pop_BZM,pop_BFE,lowavg_avgBE,lowavg_minBE,lowavg_BZM,lowavg_BFE,lowBE_avgBE,lowBE_minBE,lowBE_BZM,lowBE_BFE,lowBZM_avgBE,lowBZM_minBE,lowBZM_BZM,lowBZM_BFE,lowBFE_avgBE,lowBFE_minBE,lowBFE_BZM,lowBFE_BFE))
        else:
            f_out.write('%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n'%(all_avgBE,all_minBE,all_BZM,all_BFE,pop_avgBE,pop_minBE,pop_BZM,pop_BFE,lowavg_avgBE,lowavg_minBE,lowavg_BZM,lowavg_BFE,lowBE_avgBE,lowBE_minBE,lowBE_BZM,lowBE_BFE,lowBZM_avgBE,lowBZM_minBE,lowBZM_BZM,lowBZM_BFE,lowBFE_avgBE,lowBFE_minBE,lowBFE_BZM,lowBFE_BFE))
    f_out.close()

    #Get final csv file with selected PELE metric and final poses (with all_lowBE)
    f_in = open('%s_pelemetrics.csv'%outname,'r')
    f_out = open('%s.csv'%outname,'w')
    if glide_csv:
        header = 'molecule,smile,batch,maxtanimoto,charge,glidegscore,all_minBE,%s\n'%metric
        if not smile:
            header = header.replace('smile,','')
        if not charge:
            header = header.replace('charge,','')
        if not batch:
            header = header.replace('batch,','')
        if not maxtanimoto:
            header = header.replace('maxtanimoto,','')
    else:
        header = 'molecule,all_minBE,%s\n'%metric
    f_out.write(header)

    if not os.path.exists('%s_poses'%outname):
        os.mkdir('%s_poses'%outname)

    for i,line in enumerate(f_in):
        #Final csv
        line = line.replace('\n','').split(',')
        if i==0:
            idx_metric = line.index(metric)
            idx_minBE = line.index('all_minBE')
            if glide_csv:
                if smile: idx_smile = line.index('smile')
                if batch: idx_batch = line.index('batch')
                if maxtanimoto: idx_maxtanimoto = line.index('maxtanimoto')
                if charge: idx_charge = line.index('charge')
                idx_glidegscore = line.index('glidegscore')
            continue
        molecule = line[0]
        print('-', i, molecule)
        f_out.write('%s,'%molecule)
        if glide_csv:
            if smile: f_out.write('%s,'%line[idx_smile])
            if batch: f_out.write('%s,'%line[idx_batch])
            if maxtanimoto: f_out.write('%s,'%line[idx_maxtanimoto])
            if charge: f_out.write('%s,'%line[idx_charge])
            f_out.write('%s,'%line[idx_glidegscore])
        f_out.write('%s,%s\n'%(line[idx_minBE],line[idx_metric]))

        #Get poses
        poses_path = '%s%s/analysis/top_poses'%(pele_output,molecule)
        poses = glob.glob('%s/*.pdb'%poses_path)
        poses_BE = np.asarray([float(pose.split('BindEner')[1].split('_AtomDist')[0].split('.pdb')[0]) for pose in poses])
        #poses_BE = np.asarray([float(pose.split('BindEner')[1].split('.pdb')[0]) for pose in poses])
        idx_minBE_pose = np.argmin(poses_BE)
        cmd = 'cp %s %s_poses/%s.pdb'%(poses[idx_minBE_pose],outname,molecule)
        os.system(cmd)

        #Add connectivity
        with open('%s_poses/%s.pdb'%(outname,molecule),'r') as f:
            lines = f.readlines()
        with open('%s_poses/%s.pdb'%(outname,molecule),'w') as f:
            for line in lines:
                if 'END' not in line:
                    f.write(line)
        cmd = 'python scripts/recover_connectivity.py %s%s_prep.pdb --pdb_target %s_poses/%s.pdb -lc L -lt L'%(ligs_dir,molecule,outname,molecule)
        os.system(cmd)
    f_in.close()
    f_out.close()
