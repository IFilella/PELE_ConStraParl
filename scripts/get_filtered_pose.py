import numpy as np
import pandas as pd
import argparse
import copy
import mdtraj as mdt
import os

outnames = []

def extract_pdb_from_filter(dfaux, resdir):
    xtc_files = dfaux['trajectory']
    step = dfaux['numberOfAcceptedPeleSteps'].tolist()
    top_file = '%s/output/topologies/topology_0.pdb'%resdir    
    for i, xtc_file in enumerate(xtc_files):
        xtc_file = resdir + xtc_file
        traj = mdt.load(xtc_file, top=top_file)
        frame = int(step[i])
        struct = traj[frame]
        outname = xtc_file.replace('.xtc', '').split('/')
        outname = 'epoch_' + outname[-2] + '_' + outname[-1] + '_step_' + str(frame)
        outnames.append(outname)
        struct.save('%s/filtered_poses/%s.pdb'%(resdir, outname))
 
def main(data,fnames,thresh_up,thresh_low, resdir):
	df = pd.read_csv(data)
	dfaux = copy.copy(df)
	for i,fname in enumerate(fnames):
		print('num filter low up')
		print(i,fname,thresh_low[i],thresh_up[i])
		dfaux = dfaux[(dfaux[fname] < thresh_up[i])]
		dfaux = dfaux[(dfaux[fname] > thresh_low[i])]
		print(dfaux.shape)
	print(dfaux)
	
	# store the results: csv and pdb poses
	if not os.path.exists('%s/filtered_poses'%resdir):
            os.mkdir('%s/filtered_poses'%resdir)
	extract_pdb_from_filter(dfaux, resdir)
	if dfaux.empty == False:
		dfaux['outname'] = outnames
		csv_name = '%s_%s_%s'%(fnames, thresh_up, thresh_low)
		csv_name = csv_name.replace('[','').replace(']','').replace('\'','').replace(' ','')
		dfaux.to_csv('%s/filtered_poses/filt_data_%s.csv'%(resdir,csv_name), index=False)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filtering data.csv from results')
	parser.add_argument('-i', dest='data', help='Data csv file', required=True)
	parser.add_argument('--resdir', dest='resdir', help='Directory where the PELE results are located', required=True)
	parser.add_argument('--fnames', dest='fnames', nargs='+' , help='data.csv columns to use as filters', required = True)
	parser.add_argument('--up_t', dest='thresh_up', nargs='+', type=float, help='list containing upper thresholds (it must be as long as fnames list)')
	parser.add_argument('--low_t', dest='thresh_low', nargs='+', type=float,help='list containing lower thresholds (it must be as long as fnames list)')
	
	args = parser.parse_args()
	data = args.data
	resdir = args.resdir
	fnames = args.fnames
	thresh_up = args.thresh_up
	thresh_low = args.thresh_low
	if len(fnames) != len(thresh_up) or len(fnames) != len(thresh_low) or len(thresh_up) != len(thresh_low):
		raise ValueError('Filter names and thresholds must have the same length')
	
	main(data,fnames,thresh_up,thresh_low,resdir)
