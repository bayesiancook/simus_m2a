#! /usr/bin/python3.5

from prepare_subexperiment import *
from runcodeml import * 
from runm2a import * 
from runmm2a import * 

path2codeml_occigen = "/home/nlartillot/paml4.8/src/"
path2codeml_pbil = "/beegfs/home/lartillot/othersofts/paml4.8/src/"
path2bayescode_occigen = "/scratch/cnt0027/bbe0449/nlartillot/bayescode/data/"

exp_folder = sys.argv[1] 

# codeml

runcodeml(exp_folder, machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 24, mem=16, path2batch = exp_folder + "/codeml/", path2run = path2codeml_occigen)

# single gene bayescode: on occigen

uninfm2a_options= "-x 1 200 -pi 1.0 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi50m2a_options= "-x 1 200 -pi 0.5 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"

runm2a(exp_folder, uninfm2a_options, "uninfm2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = exp_folder + "/singlegene/", path2run = path2bayescode_occigen)

runm2a(exp_folder, uninfpi50m2a_options, "uninfpi50m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = exp_folder + "/singlegene/", path2run = path2bayescode_occigen)

# multi gene bayescode: on occigen

indmm2a_options= "-x 1 200 -nucrates ind -bl ind +G"
shrunkenmm2a_options= "-x 1 200 -nucrates shrunken -bl shrunken +G"

runmm2a(exp_folder, indmm2a_options, "indmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = exp_folder + "/multigene/", path2run = path2bayescode_occigen)
runmm2a(exp_folder, shrunkenmm2a_options, "shrunkenmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = exp_folder + "/multigene/", path2run = path2bayescode_occigen)

