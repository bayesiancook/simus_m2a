#! /usr/bin/python3.5

import postpredsimu
import sys

# original subexperiment on which to do the post pred simulation should 
exp_folder = sys.argv[1]
# simulation based on uninfm2a chains in singlegene folder
basename = "uninfm2a"

# list of genes in decreasing order of delta lnL (codeml)
listname = sys.argv[2]

# target proportion of genes simulated under positive selection
prop_pos = float(sys.argv[3])

# shrinking posw and dposom
shrink_posw = float(sys.argv[4])
shrink_dposom = float(sys.argv[5])

# name of target directory
target_name = sys.argv[6]

postpredsimu.single_gene_postpred(exp_folder, listname, basename, target_name, prop_pos = prop_pos, shrink_posw = shrink_posw, shrink_dposom = shrink_dposom, prob_posw = 1.0, prob_dposom = 1.0, burnin=100)

