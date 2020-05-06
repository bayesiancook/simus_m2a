#! /usr/bin/python3.5

#
# post analysis of experiments on original empirical data
# 

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
emp_folder = exp_folder + "/empirical"
res_folder = exp_folder + "/results"

if not os.path.exists(res_folder):
    print("creating results folder")
    os.system("mkdir " + res_folder)

# post analysis

single_basename = ["uninfm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a", "unconsindmm2a", "unconsshrunkenmm2a"]

m2a_postanalysis(emp_folder, single_basename, multi_basename, outname = exp_folder + "/results/emp")

