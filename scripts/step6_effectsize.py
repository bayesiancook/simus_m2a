#! /usr/bin/python3.5

# a restricted post analysis under codeml and multi gene models

import sys
import os

from postes import * 

exp_folder = sys.argv[1]

simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_posw03", "simushrunken"]

multi_basename = ["unconsindmm2a", "indmm2a", "unconsshrunkenmm2a", "shrunkenmm2a"]

full_m2a_postes(exp_folder, ["simu10"], ["uninfm2a"], multi_basename, "allsimupos", esmode = "pos", cutoff_list = [0.1, 0.2, 0.3])

# full_m2a_postes(exp_folder, simu_list, [], multi_basename, "allsimupos", esmode = "pos", cutoff_list = [0.1, 0.2, 0.3])




# full_m2a_postes(exp_folder, simu_list, [], multi_basename, "allsimuexcess", esmode = "excess", cutoff_list = [0.05, 0.1, 0.2])

