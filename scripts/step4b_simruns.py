#! /usr/bin/python3.5

import sys
import os

from simuparams import * 
from runmm2a import * 

path2bayescode_p2chpd = "/home_nfs/lbbe/nlartillot/bayescode/data/"

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"

simu_list = ["simu30", "simu10", "simu03"]
# simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01", "simushared", "simushrunken"]

# prepare runs
# for all simus: multigene oracle + uninf 
# for all simus: codeml and shared shrunken ind multigene
# for all single-gene simulations: also uncons and uninf multigene
# for simu30 10 and 03: single gene m2a analyses

gene_file = exp_dir + "/empirical/all.list"
with open(gene_file, 'r') as infile:
    ngene = len([line for line in infile])

options = {
        "uninfpi50mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.50 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "uninfpi10mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.10 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "uninfpi02mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.02 -purom 0.5 0.5 -dposom 10 1 -purw 0.5 0.5 -posw 0.5 0.5",
        "subjpi50mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.50 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1",
        "subjpi10mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.10 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1",
        "subjpi02mm2a" : "-x 1 1100 -nucrates ind -bl ind +G -pi 0.02 -purom 0.5 0.5 -dposom 2 1 -purw 0.5 0.5 -posw 0.1 0.1"
}

for simu in simu_list:

    simu_folder = exp_dir + simu 

    multi_dir = simu_folder + "/multigene2/"
    if not os.path.exists(simu_folder + "/multigene2"):
        os.system("mkdir " + simu_folder + "/multigene2")
        os.system("cp " + simu_folder + "/multigene/*ali " + simu_folder + "/multigene2/")
        os.system("cp " + simu_folder + "/multigene/all.* " + simu_folder + "/multigene2/")

    # multi gene bayescode
    core = 16
    nodes = 3

    for name in options:
        runmm2a(simu_folder, options[name], name, machine="p2chpd", queue="parallel", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene2/", path2run = path2bayescode_p2chpd)

    (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = get_empirical_moments(simu_folder)

    infmm2a_options= "-x 1 1100 -nucrates ind -bl ind +G -pi {0} -purom {1} {2} -dposom {3} {4} -purw {5} {6} -posw {7} {8}".format(pi, purom_mean, purom_invconc, dposom_mean, dposom_invshape, purw_mean, purw_invconc, posw_mean, posw_invconc)

    runmm2a(simu_folder, infmm2a_options, "infmm2a", machine="p2chpd", queue="parallel", nodes=nodes, core=core, time=24, mem=16, path2batch = simu_folder + "/multigene2/", path2run = path2bayescode_p2chpd)


