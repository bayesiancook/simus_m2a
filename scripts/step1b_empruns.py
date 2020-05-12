#! /usr/bin/python3.5

#
# experiments on original empirical data
# also: prepare ppred scripts for multigene data
# 

import sys
import os

from new_experiment import * 
from prepare_subexperiment import *
from runcodeml import * 
from runm2a import * 
from runmm2a import * 

# input:
data_folder = sys.argv[1]
gene_file = sys.argv[2]
tree_file = sys.argv[3]
exp_folder = sys.argv[4]
taxon_file = sys.argv[5]

# machines on which to run non-parallel jobs
machine = "pbil"
queue = "none"
path2codeml = "/beegfs/home/lartillot/othersofts/paml4.8/src/"
path2bayescode = "/beegfs/home/lartillot/code/bayescode/data/"

# machines on which to run parallel jobs
multi_machine = "p2chpd"
multi_queue = "parallel"
multi_path2bayescode = "/home_nfs/lbbe/nlartillot/bayescode/data/"

# multi_machine = "occigen"
# multi_queue = "none"
# multi_path2bayescode = "/scratch/cnt0027/bbe0449/nlartillot/bayescode/data/"

# make new experiment 
new_experiment(data_folder, gene_file, taxon_file, tree_file, exp_folder)

# get gene list
with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

# determine number of nodes and cores for parallel runs
# for p2chpd
core = 16
nodes = 3

emp_folder = exp_folder + "/empirical"

# set up subexp on original_data
prepare_subexperiment(emp_folder)

# prepare all runs

name2command = {"sharedmm2a" : "-x 1 1100 -nucrates shared -bl shared +G",
        "shrunkenmm2a" : "-x 1 1100 -nucrates shrunken -bl shrunken +G",
        "indmm2a" : "-x 1 1100 -nucrates ind -bl ind +G"}

for name in name2command:
    runmm2a(emp_folder, name2command[name], name, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)

# multigene post pred bash scripts (prepare them already at that step)
for name in ["sharedmm2a", "shrunkenmm2a", "indmm2a"]:
    runppredmm2a(emp_folder, name, burnin = 500, every = 30, nrep = 1, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=1, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)


