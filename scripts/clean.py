import os
import glob

dirs = glob.glob('R0_test/*/')
for d in dirs:
    cmd = 'rm -r %soutput/'%d
    print(cmd)
    os.system(cmd)
    cmd = 'rm -r %sanalysis_1/'%d
    print(cmd)
    os.system(cmd)
    cmd = 'rm -r %sresults/top_poses/'%d
    os.system(cmd)
    cmd = 'rm -r %sresults_1/'%d
    os.system(cmd)
