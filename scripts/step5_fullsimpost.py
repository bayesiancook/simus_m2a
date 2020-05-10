#! /usr/bin/python3.5

# a restricted post analysis under codeml and multi gene models

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
# fields = ["n", "fdr"]
fields = ["n", "fdr", "efnr", "fnr", "ees", "es", "e%es", "%es"]

simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01", "simushrunken", "simushared"]

multi_basename = ["unconsindmm2a", "indmm2a", "unconsshrunkenmm2a", "shrunkenmm2a", "sharedmm2a"]

full_m2a_postanalysis(exp_folder, ["simu30", "simu10", "simu30_shrink_dposom03", "simu30_shrink_posw03"], [], ["unconsindmm2a"], "allsimu", fields = fields, with_tex = True)
# full_m2a_postanalysis(exp_folder, simu_list, [], multi_basename, "allsimu", fields = fields, with_tex = True)

