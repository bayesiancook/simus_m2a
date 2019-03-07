#! /usr/bin/python3.5

import postanalysis
import sys

exp_folder = sys.argv[1]
outname = sys.argv[2]

# single_basename = ["uninfm2a", "uninfpi50m2a", "uninfpi20m2a", "uninfpi02m2a", "infpi50m2a", "infpi20m2a", "infm2a", "infpi02m2a"]
# multi_basename = ["infmm2a", "indmm2a", "shrunkenmm2a"]
single_basename = ["uninfm2a", "uninfpi50m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a"]

postanalysis.m2a_postanalysis(exp_folder, single_basename, multi_basename, outname = outname, burnin = 20, dlnlmin = 3, min_omega = 1.2)

