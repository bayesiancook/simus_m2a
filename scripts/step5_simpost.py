#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
res_dir = exp_dir + "results/"

# post analysis 
single_basename = ["infm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a"]
full_multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a", "unconsindmm2a", "unconsshrunkenmm2a", "uninfindmm2a", "uninfshrunkenmm2a"]

print("single gene analyses")
for simu in ["simu30", "simu10", "simu03"]:
    print(simu)
    m2a_postanalysis(exp_dir + simu, single_basename, [], outname = res_dir + "single" + simu, burnin = 100)

print("extended multi gene analyses")
for simu in ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01"]:
    print(simu)
    m2a_postanalysis(exp_dir + simu, [], full_multi_basename, outname = res_dir + "multi" + simu, burnin = 100)

print("restricted multi gene analyses")
for simu in ["simushared", "simushrunken"]:
    print(simu)
    m2a_postanalysis(exp_dir + simu, [], multi_basename, outname = res_dir + simu, burnin = 100)

