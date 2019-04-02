#! /usr/bin/python3.5

import sys
import os

from simuparams import * 
from runcodeml import * 
from runm2a import * 
from runmm2a import * 

path2codeml_occigen = "/home/nlartillot/paml4.8/src/"
path2codeml_pbil = "/beegfs/home/lartillot/othersofts/paml4.8/src/"
path2bayescode_occigen = "/scratch/cnt0027/bbe0449/nlartillot/bayescode/data/"

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"

simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01", "simushared", "simushrunken"]

# prepare runs
# for all simus: codeml and shared shrunken ind multigene
# for all single-gene simulations: also uncons and uninf multigene
# for simu30 10 and 03: single gene m2a analyses

gene_file = exp_dir + "/empirical/all.list"
with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

for simu in simu_list:

    simu_folder = exp_dir + simu 

    # codeml
    runcodeml(exp_dir + simu, machine="pbil", queue="none", core = 1, njobs_per_batch = 1, time = 24, mem=2, path2batch = simu_folder + "/codeml/", path2run = path2codeml_pbil)
    # runcodeml(exp_dir + simu, machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 24, mem=16, path2batch = simu_folder + "/codeml/", path2run = path2codeml_occigen)

    # multi gene bayescode

    core = 24
    nodes = ngene // core // 5

    sharedmm2a_options= "-x 1 1100 -nucrates shared -bl shared +G"
    shrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G"
    indmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G"

    unconsshrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G -unconsprior"
    unconsindmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G -unconsprior"

    uninfshrunkenmm2a_options= "-x 1 1100 -pi 0.5 0.5 -nucrates shrunken -bl shrunken +G -unconsprior"
    uninfindmm2a_options= "-x 1 1100 -pi 0.5 0.5 -nucrates ind -bl ind +G -unconsprior"

    runmm2a(simu_folder, sharedmm2a_options, "sharedmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, shrunkenmm2a_options, "shrunkenmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, indmm2a_options, "indmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)

    runmm2a(simu_folder, unconsshrunkenmm2a_options, "unconsshrunkenmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, unconsindmm2a_options, "unconsindmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)

    runmm2a(simu_folder, uninfshrunkenmm2a_options, "uninfshrunkenmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, uninfindmm2a_options, "uninfindmm2a", machine="occigen", queue="none", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)


