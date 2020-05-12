#! /usr/bin/python3.5

# more generic than runm2a (does not assume the structure of a simulation experiment)
# assumes single genes have each their own tree

import sys
import os
from ete3 import Tree
from datalib import *
from makebatch import *

path2bayescode = "/beegfs/home/lartillot/code/bayescode/data/"

def runsinglegene(exp_folder, gene_list_name, options, basename, batch_mode = "sbatch", machine = "pbil", queue = "none", core = 1, njobs_per_batch = 1, time = 12, mem=1, path2run = path2bayescode):

    exp_dir = exp_folder + "/"
    if not os.path.exists(exp_dir):
        print("directory does not exist")
        sys.exit()

    # get gene list
    with open(exp_dir + gene_list_name) as file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in file]

    # make batch file 
    command_list = ["codonm2a -d " + gene + ".ali -t " + gene + ".tree " + options + " " + basename + gene for gene in gene_list]
    makebatchlist(command_list, basename, time=time, mem=mem, machine=machine, queue=queue, core=core, njobs_per_batch=njobs_per_batch, mode = batch_mode, path2batch = exp_dir, path2run = path2run)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("runsinglegene.py exp_folder gene_list options basename")
        sys.exit()

    runsinglegene(*sys.argv[1:])


