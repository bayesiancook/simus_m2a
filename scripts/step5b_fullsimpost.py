#! /usr/bin/python3.5

# a restricted post analysis under codeml and multi gene models

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
fields = ["n", "fdr", "efnr", "fnr"]

simu_list = ["simushrunken", "simushared", "simuind"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a"]
full_m2a_postanalysis(exp_folder, simu_list, [], multi_basename, "allsimu", fields = fields)

