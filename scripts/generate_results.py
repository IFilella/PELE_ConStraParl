from pele_platform.analysis import Analysis
import argparse
import glob
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='')
    parser.add_argument('--workDIR', dest="workDIR", help = "",required=True)
    args = parser.parse_args()

    #Parse inputs
    workDIR = args.workDIR

    originalDIR = os.getcwd()

    os.chdir(workDIR)
    path = '/'.join(workDIR.split('/')[:-2])
    path = os.path.join(path,'2_pele_pdb_preprocessor/complexes')
    topology = path + '/complex_2.pdb'
    analysis = Analysis(resname="LIG", chain="L", simulation_output="output",traj='trajectory.xtc',topology=topology,be_column=4,report="report",cpus=32)
    analysis.generate(path="analysis", clustering_type="meanshift")
    os.chdir(originalDIR)
