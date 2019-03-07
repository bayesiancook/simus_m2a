#! /usr/bin/python3.5

import sys
import os
from ete3 import Tree
from datalib import *
from makebatch import *

def runm2a(exp_folder, options, basename, batch_mode = "sbatch", machine = "p2chpd", queue = "parallel2", core = 8, njobs_per_batch = 8, time = 50, mem=16, path2batch = "", path2run = ""):

    exp_dir = exp_folder + "/"
    data_dir = exp_dir + "data/"
    single_dir = exp_dir + "singlegene/"
    if not os.path.exists(single_dir):
        print("singlegene directory does not exist")
        sys.exit()

    # get gene list
    with open(exp_dir + "all.list") as file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in file]

    # make batch file 
    command_list = ["codonm2a -d " + gene + ".ali -t " + gene + ".tree " + options + " " + basename + gene for gene in gene_list]
    makebatchlist(command_list, basename, time=time, mem=mem, machine=machine, queue=queue, core=core, njobs_per_batch=njobs_per_batch, mode = batch_mode, path2batch = path2batch, path2run = path2run)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("runm2a.py experiment/subexperiment options basename")
        sys.exit()

    runm2a(*sys.argv[1:])


