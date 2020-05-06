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

single_simu_list = ["simu30", "simu10", "simu03"]
simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01", "simushared", "simushrunken"]

# prepare runs
# for all simus: codeml and shared shrunken ind unconsshrunken unconsind multigene
# for simu30 10 and 03: single gene m2a analyses

gene_file = exp_dir + "/empirical/all.list"
with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

# for p2chpd
core = 16
nodes = 3

# for occigen
# core = 24
# nodes = ngene // core // 5

# single-gene bayescode options
# uninfm2a_options= "-x 1 600 -pi 1.0 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
uninfpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
subjpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
subjpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
subjpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"

# multigene bayescode options
sharedmm2a_options= "-x 1 1100 -nucrates shared -bl shared +G"
shrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G"
indmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G"
unconsshrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G -unconsprior"
unconsindmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G -unconsprior"

for simu in single_simu_list:

    simu_folder = exp_dir + simu 

    runm2a(simu_folder, uninfm2a_options, "uninfm2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, uninfpi50m2a_options, "uninfpi50m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, uninfpi10m2a_options, "uninfpi10m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, uninfpi02m2a_options, "uninfpi02m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, subjpi50m2a_options, "subjpi50m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, subjpi10m2a_options, "subjpi10m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)
    runm2a(simu_folder, subjpi02m2a_options, "subjpi02m2a", machine=machine, queue=queue, core=1, njobs_per_batch=1, time = 6, mem=1, path2batch = emp_folder + "/singlegene/", path2run=path2bayescode)

for simu in simu_list:

    simu_folder = exp_dir + simu 

    # codeml
    runcodeml(exp_dir + simu, machine="pbil", queue="none", core = 1, njobs_per_batch = 1, time = 24, mem=2, path2batch = simu_folder + "/codeml/", path2run = path2codeml_pbil)

    runmm2a(simu_folder, sharedmm2a_options, "sharedmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = multi_path2bayescode)
    runmm2a(simu_folder, shrunkenmm2a_options, "shrunkenmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = multi_path2bayescode)
    runmm2a(simu_folder, indmm2a_options, "indmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = multi_path2bayescode)
    runmm2a(simu_folder, unconsshrunkenmm2a_options, "unconsshrunkenmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = multi_path2bayescode)
    runmm2a(simu_folder, unconsindmm2a_options, "unconsindmm2a", machine=multi_machine, queue=multi_queue, nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = multi_path2bayescode)

