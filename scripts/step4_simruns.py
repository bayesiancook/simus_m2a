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
single_simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01"]
multi_simu_list = ["simushared", "simushrunken"]

# prepare runs
# for all simus: codeml and shared shrunken ind multigene
# for all single-gene simulations: also uncons and uninf multigene
# for simu30 10 and 03: single gene m2a analyses

for simu in simu_list:

    # codeml
    runcodeml(exp_dir + simu, machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 24, mem=16, path2batch = simu_folder + "/codeml/", path2run = path2codeml_occigen)

    # multi gene bayescode

    sharedmm2a_options= "-x 1 1100 -nucrates shared -bl shared +G"
    shrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G"
    indmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G"

    runmm2a(simu_folder, sharedmm2a_options, "sharedmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, shrunkenmm2a_options, "shrunkenmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, indmm2a_options, "indmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)


for simu in single_simu_list:

    # additional multi gene bayescode runs

    unconsshrunkenmm2a_options= "-x 1 1100 -nucrates shrunken -bl shrunken +G -unconsprior"
    unconsindmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G -unconsprior"
    uninfshrunkenmm2a_options= "-x 1 1100 -pi 0.5 0.5 -nucrates shrunken -bl shrunken +G -unconsprior"
    uninfindmm2a_options= "-x 1 1100 -pi 0.5 0.5 -nucrates ind -bl ind +G -unconsprior"

    runmm2a(simu_folder, unconsshrunkenmm2a_options, "unconsshrunkenmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, unconsindmm2a_options, "unconsindmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)

    runmm2a(simu_folder, uninfshrunkenmm2a_options, "uninfshrunkenmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)
    runmm2a(simu_folder, uninfindmm2a_options, "uninfindmm2a", machine="occigen", queue="none", nodes=2, core=24, time=12, mem=16, path2batch = simu_folder + "/multigene/", path2run = path2bayescode_occigen)

# single gene bayescode runs

for simu in ["simu30", "simu10", "simu03"]:

    simu_folder = exp_folder + "/simu30"

    uninfpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
    uninfpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"
    uninfpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5"

    subjpi50m2a_options= "-x 1 600 -pi 0.50 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
    subjpi10m2a_options= "-x 1 600 -pi 0.10 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
    subjpi02m2a_options= "-x 1 600 -pi 0.02 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"

    (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = get_empirical_moments(simu_folder)

    infm2a_options= "-x 1 600 -pi {0} -purom {1} {2} -dposom {3} {4} -purw {5} {6} -posw {7} {8}".format(pi, purom_mean, purom_invconc, dposom_mean, dposom_invshape, purw_mean, purw_invconc, posw_mean, posw_invconc)

    runm2a(simu_folder, uninfpi50m2a_options, "uninfpi50m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)
    runm2a(simu_folder, uninfpi10m2a_options, "uninfpi10m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)
    runm2a(simu_folder, uninfpi02m2a_options, "uninfpi02m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)

    runm2a(simu_folder, subjpi50m2a_options, "subjpi50m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)
    runm2a(simu_folder, subjpi10m2a_options, "subjpi10m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)
    runm2a(simu_folder, subjpi02m2a_options, "subjpi02m2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)

    runm2a(simu_folder, infm2a_options, "infm2a", machine="occigen", queue="none", core = 24, njobs_per_batch = 24, time = 2, mem=16, path2batch = simu_folder + "/singlegene/", path2run = path2bayescode_occigen)

