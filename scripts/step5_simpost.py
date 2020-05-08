#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
res_dir = exp_dir + "results/"
with_sites = False

# an extended post analysis for all methods (codeml, single gene and multigene) for simu30 / 10 / 03

single_basename = ["infm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]

multi_basename = ["unconsindmm2a", "unconsshrunkenmm2a", "indmm2a", "shrunkenmm2a", "sharedmm2a"]

simu_list = ["simu10"]
# simu_list = ["simu30", "simu10", "simu03"]

for simu in simu_list:
    print(simu)
    m2a_postanalysis(exp_dir + simu, single_basename, multi_basename, outname = res_dir + "single" + simu, with_sites = with_sites)

