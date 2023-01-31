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
    topology = glob.glob('%s/input/*processed*'%workDIR)[0]
    analysis = Analysis(resname="LIG", chain="L", simulation_output="output",traj='trajectory.xtc',topology=topology,be_column=4,report="report",cpus=32)
    analysis.generate(path="results", clustering_type="meanshift")
    os.chdir(originalDIR)
