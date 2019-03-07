#! /usr/bin/python3.5

import sys
import os
from ete3 import Tree
from datalib import *
from makebatch import *

def runmm2a(exp_folder, options, basename, machine = "p2chpd", queue = "parallel2", nodes=1, core = 8, time = 50, mem=16, batch_mode = "sbatch", path2batch = "", path2run = ""):

    exp_dir = exp_folder + "/"
    data_dir = exp_dir + "data/"
    multi_dir = exp_dir + "multigene/"
    if not os.path.exists(multi_dir):
        print("multigene directory does not exist")
        sys.exit()

    # make batch file
    command = "multigenecodonm2a -d all.list -t all.tree " + options + " " + basename 
    makebatch(command, basename, time=time, mem=mem, machine=machine, queue=queue, nodes=nodes, core=core, mode = batch_mode, path2batch = path2batch, path2run=path2run)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("runmm2a.py experiment/subexperiment options basename")
        sys.exit()

    runmm2a(*sys.argv[1:])


