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

taxon_file = "alltaxa"

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
# for occigen
# core = 24
# nodes = ngene // core // 5

emp_folder = exp_folder + "/empirical"

# set up subexp on original_data
prepare_subexperiment(emp_folder)

# prepare all runs

# codeml
# make scripts for job submission
runcodeml(emp_folder, machine=machine, queue=queue, core=1, njobs_per_batch=1, time=24, mem=1, path2batch=emp_folder + "/codeml/", path2run=path2codeml)

# single gene bayescode 

name2command = {"uninfm2a" : "-x 1 600 -pi 1.0 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "uninfpi50m2" : "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "uninfpi10m2" : "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "uninfpi02m2" : "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "subjpi50m2a" : "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1",
        "subjpi10m2a" : "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1",
        "subjpi02m2a" : "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"}

for name in name2command:
    runm2a(emp_folder, name2command[name], name, machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)

# multi gene bayescode

name2command = {"sharedmm2a" : "-x 1 1100 -nucrates shared -bl shared +G",
        "shrunkenmm2a" : "-x 1 1100 -nucrates shrunken -bl shrunken +G",
        "indmm2a" : "-x 1 1100 -nucrates ind -bl ind +G",
        "unconsshrunkenmm2a" : "-x 1 1100 -nucrates shrunken -bl shrunken +G -unconsprior",
        "unconsindmm2a" : "-x 1 1100 -nucrates ind -bl ind +G -unconsprior"}

for name in name2command:
    runmm2a(emp_folder, name2command[name], name, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)

# multigene post pred bash scripts (prepare them already at that step)
for name in ["sharedmm2a", "shrunkenmm2a"]:
    runppredmm2a(emp_folder, name, burnin = 500, every = 30, nrep = 1, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=1, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)


