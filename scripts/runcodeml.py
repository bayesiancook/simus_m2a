#! /usr/bin/python3.5

import sys
import os
from ete3 import Tree
from datalib import *
from makebatch import *

path2ctl = "/beegfs/home/lartillot/othersofts/paml4.8/"

def runcodeml(exp_folder, machine="p2chpd", batch_mode = "sbatch", queue = "parallel2", core = 8, njobs_per_batch = 8, time = 50, mem=16, path2batch = "", path2run = ""):

    exp_dir = exp_folder + "/"
    codeml_dir = exp_dir + "codeml/"
    if not os.path.exists(codeml_dir):
        print("directory " + codeml_dir + " does not exist")
        sys.exit()

    # get codeml control file
    with open(codeml_dir + "yes", 'w') as stdfile:
        stdfile.write("\n")
        stdfile.write("yes\n")
        stdfile.write("\n")

    with open(path2ctl + "ctl0", 'r') as ctlfile:
        ctl0 = ctlfile.read()

    # get gene list
    with open(exp_dir + "all.list") as file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in file]

    # get tree
    tree = Tree(exp_dir + "all.tree")

    # prepare control files (with .ctl extension) for all genes
    for gene in gene_list:

        with open(codeml_dir + gene + ".ctl", 'w') as ctlfile:
            ctlfile.write("       seqfile =  " + gene + ".ali\n")
            ctlfile.write("      treefile =  " + gene + ".tree\n")
            ctlfile.write("       outfile =  " + gene + ".codeml\n")
            ctlfile.write(ctl0)

    # make batch file 
    command_list = ["codeml " + gene + ".ctl < yes" for gene in gene_list]
    makebatchlist(command_list, "codeml_analysis", machine=machine, time=time, queue=queue, mem=mem, core=core, njobs_per_batch=njobs_per_batch, mode = batch_mode, path2batch=path2batch, path2run=path2run)

def reruncodeml(exp_folder, machine="p2chpd", batch_mode = "sbatch", queue = "parallel2", core = 8, njobs_per_batch = 8, time = 50, mem=16, path2batch = "", path2run = ""):

    exp_dir = exp_folder + "/"
    codeml_dir = exp_dir + "codeml/"
    if not os.path.exists(codeml_dir):
        print("directory " + codeml_dir + " does not exist")
        sys.exit()

    current_dir = os.getcwd()
    os.chdir(codeml_dir)
    os.system("grep -c \"Time used\" *codeml | grep codeml:0 > rerun.list")
    with open("rerun.list", 'r') as listfile:
        gene_list = [gene.rstrip('\n').replace(".codeml:0","") for gene in listfile]

    print("rerunning:")
    for gene in gene_list:
        print(gene)
        os.system("rm " + gene + ".codeml")

    os.chdir(current_dir)

    # make batch file 
    command_list = ["codeml " + gene + ".ctl < yes" for gene in gene_list]
    makebatchlist(command_list, "rerun_codeml_analysis", machine=machine, time=time, queue=queue, mem=mem, core=core, njobs_per_batch=njobs_per_batch, mode = batch_mode, path2batch=path2batch, path2run=path2run)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("runcodeml.py experiment/subexperiment")
        sys.exit()

    runcodeml(*sys.argv[1:])

