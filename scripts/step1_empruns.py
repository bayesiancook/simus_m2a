#! /usr/bin/python3.5

import sys
import os

from new_experiment import * 
from prepare_subexperiment import *
from runcodeml import * 
from runm2a import * 
from runmm2a import * 

path2codeml_occigen = "/home/nlartillot/paml4.8/src/"
path2codeml_pbil = "/beegfs/home/lartillot/othersofts/paml4.8/src/"
path2bayescode_occigen = "/scratch/cnt0027/bbe0449/nlartillot/bayescode/data/"

# input:

data_folder = sys.argv[1]
gene_file = sys.argv[2]
tree_file = sys.argv[3]
exp_folder = sys.argv[4]

taxon_file = "alltaxa"

# make new experiment 
new_experiment(data_folder, gene_file, taxon_file, tree_file, exp_folder)

with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])
#
# sub-experiments on original empirical data
# 

emp_folder = exp_folder + "/empirical"

# set up subexp on original_data
prepare_subexperiment(emp_folder)

# prepare all runs (on occigen)

# codeml
runcodeml(emp_folder, machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 24, mem=16, path2batch = emp_folder + "/codeml/", path2run = path2codeml_occigen)

# single gene bayescode 
uninfm2a_options= "-x 1 600 -pi 1.0 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
subjpi10m2a_options= "-x 1 600 -pi 0.1 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"

runm2a(emp_folder, uninfm2a_options, "uninfm2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 6, mem=16, path2batch = emp_folder + "/singlegene/", path2run = path2bayescode_occigen)
runm2a(emp_folder, uninfpi50m2a_options, "uninfpi50m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 6, mem=16, path2batch = emp_folder + "/singlegene/", path2run = path2bayescode_occigen)
runm2a(emp_folder, uninfpi10m2a_options, "uninfpi10m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 6, mem=16, path2batch = emp_folder + "/singlegene/", path2run = path2bayescode_occigen)
runm2a(emp_folder, uninfpi02m2a_options, "uninfpi02m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 6, mem=16, path2batch = emp_folder + "/singlegene/", path2run = path2bayescode_occigen)
runm2a(emp_folder, subjpi10m2a_options, "subjpi10m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 6, mem=16, path2batch = emp_folder + "/singlegene/", path2run = path2bayescode_occigen)

# multi gene bayescode
sharedmm2a_options= "-x 1 1100 -nucrates shared -bl shared +G"
shrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G"
indmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G"

core = 24
nodes = ngene // core // 5

runmm2a(emp_folder, sharedmm2a_options, "sharedmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = emp_folder + "/multigene/", path2run = path2bayescode_occigen)
runmm2a(emp_folder, shrunkenmm2a_options, "shrunkenmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = emp_folder + "/multigene/", path2run = path2bayescode_occigen)
runmm2a(emp_folder, indmm2a_options, "indmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = emp_folder + "/multigene/", path2run = path2bayescode_occigen)

# multigene post pred bash scripts (prepare them already at that step)

runppredmm2a(emp_folder, "sharedmm2a", burnin = 100, every = 30, nrep = 1, machine="occigen", queue="none", nodes=nodes, core=core, time=1, mem=16, path2batch = emp_folder + "/multigene/", path2run = path2bayescode_occigen)
runppredmm2a(emp_folder, "shrunkenmm2a", burnin = 100, every = 30, nrep = 1, machine="occigen", queue="none", nodes=nodes, core=core, time=1, mem=16, path2batch = emp_folder + "/multigene/", path2run = path2bayescode_occigen)


