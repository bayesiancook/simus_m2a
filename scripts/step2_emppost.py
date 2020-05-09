#! /usr/bin/python3.5

#
# post analysis of experiments on original empirical data
# 

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]

# fields = ["ndisc"]
fields = ["ndisc", "e-fnr"]

if not os.path.exists(exp_folder + "/results/"):
    print("creating results folder")
    os.system("mkdir " + exp_folder + "/results")

# post analysis

single_basename = ["uninfm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a", "unconsindmm2a", "unconsshrunkenmm2a"]

full_m2a_postanalysis(exp_folder, ["empirical"], single_basename, multi_basename, "emp", fields = fields)

