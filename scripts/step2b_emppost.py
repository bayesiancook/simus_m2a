#! /usr/bin/python3.5

#
# post analysis of experiments on original empirical data
# 

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]

fields = ["n", "efnr"]

if not os.path.exists(exp_folder + "/results/"):
    print("creating results folder")
    os.system("mkdir " + exp_folder + "/results")

# post analysis

multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a"]

full_m2a_postanalysis(exp_folder, ["empirical"], [], multi_basename, "emp", fields = fields)

