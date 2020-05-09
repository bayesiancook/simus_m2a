#! /usr/bin/python3.5

# an extended post analysis for all methods (codeml, single gene and multigene) for simu30 / 10 / 03

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]

fields = ["ndisc", "fdr"]
# fields = ["ndisc", "fdr", "e-fnr", "fnr"]

simu_list = ["simu10"]
# simu_list = ["simu30", "simu10", "simu03"]

single_basename = ["infm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]

multi_basename = ["unconsindmm2a", "unconsshrunkenmm2a", "indmm2a", "shrunkenmm2a", "sharedmm2a"]

full_m2a_postanalysis(exp_folder, simu_list, single_basename, multi_basename, "focalsimu", fields = fields)

