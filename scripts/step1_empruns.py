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

with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

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

# options for bayescode
uninfm2a_options= "-x 1 600 -pi 1.0 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
subjpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
subjpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
subjpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"

# make scripts for job submission
runm2a(emp_folder, uninfm2a_options, "uninfm2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, uninfpi50m2a_options, "uninfpi50m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, uninfpi10m2a_options, "uninfpi10m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, uninfpi02m2a_options, "uninfpi02m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, subjpi50m2a_options, "subjpi50m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, subjpi10m2a_options, "subjpi10m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
runm2a(emp_folder, subjpi02m2a_options, "subjpi02m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)

# multi gene bayescode

# options for bayescode
sharedmm2a_options= "-x 1 1100 -nucrates shared -bl shared +G"
shrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G"
indmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G"
unconsshrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G -unconsprior"
unconsindmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G -unconsprior"

# make scripts for job submission
runmm2a(emp_folder, sharedmm2a_options, "sharedmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)
runmm2a(emp_folder, shrunkenmm2a_options, "shrunkenmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)
runmm2a(emp_folder, indmm2a_options, "indmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)
runmm2a(emp_folder, unconsshrunkenmm2a_options, "unconsshrunkenmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)
runmm2a(emp_folder, unconsindmm2a_options, "unconsindmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)

# multigene post pred bash scripts (prepare them already at that step)
runppredmm2a(emp_folder, "sharedmm2a", burnin = 500, every = 30, nrep = 1, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=1, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)
runppredmm2a(emp_folder, "shrunkenmm2a", burnin = 500, every = 30, nrep = 1, machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=1, mem=16, path2batch=emp_folder + "/multigene/", path2run=multi_path2bayescode)


