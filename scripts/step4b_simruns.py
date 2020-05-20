#! /usr/bin/python3.5

#
# codeml, single- and multi-gene bayesian analyses on simulated data
# 

import sys
import os

from simuparams import * 
from runcodeml import * 
from runm2a import * 
from runmm2a import * 

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

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"

simu_list = ["simushared", "simushrunken", "simuind"]

gene_file = exp_dir + "/empirical/all.list"
with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

# determining number of nodes and cores per node
# for p2chpd
core = 16
nodes = 3

# codeml and multi-gene runs

name2command = {"sharedmm2a" : "-x 1 1100 -nucrates shared -bl shared +G",
        "shrunkenmm2a" : "-x 1 1100 -nucrates shrunken -bl shrunken +G",
        "indmm2a" : "-x 1 1100 -nucrates ind -bl ind +G"}

for simu in simu_list:

    simu_folder = exp_dir + simu 

    for name in name2command:
        runmm2a(simu_folder, name2command[name], name, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=simu_folder + "/multigene/", path2run=multi_path2bayescode)


