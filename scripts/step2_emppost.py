#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
emp_folder = exp_folder + "/empirical"

# post analysis

single_basename = ["uninfm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a", "unconsindmm2a", "unconsshrunkenmm2a", "uninfindmm2a", "uninfshrunkenmm2a"]

m2a_postanalysis(emp_folder, single_basename, multi_basename, outname = exp_folder + "/results/emp", burnin = 100)

