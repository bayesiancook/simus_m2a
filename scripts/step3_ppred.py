#! /usr/bin/python3.5

import sys
import os

from postpredsimu import *
from prepare_subexperiment import *

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
emp_folder = exp_dir + "empirical"
postlist_name = exp_dir + "results/emp.sortedparams"

# post pred simulations, single gene
# varying proportion of genes under positive selection, and effect size (both dposom and posw)

single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu30", prop_pos = 0.30, burnin=100, shrink_posw = 1.0, shrink_dposom = 1.0, prob_posw = 1.0, prob_dposom = 1.0)
single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu10", prop_pos = 0.10, burnin=100, shrink_posw = 1.0, shrink_dposom = 1.0, prob_posw = 1.0, prob_dposom = 1.0)
single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu03", prop_pos = 0.01, burnin=100, shrink_posw = 1.0, shrink_dposom = 1.0, prob_posw = 1.0, prob_dposom = 1.0)

single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu30_shrink_dposom03", prop_pos = 0.30, burnin=100, shrink_posw = 1.0, shrink_dposom = 0.3, prob_posw = 1.0, prob_dposom = 1.0)
single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu30_shrink_dposom01", prop_pos = 0.30, burnin=100, shrink_posw = 1.0, shrink_dposom = 0.1, prob_posw = 1.0, prob_dposom = 1.0)
single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu30_shrink_posw03", prop_pos = 0.30, burnin=100, shrink_posw = 0.3, shrink_dposom = 1.0, prob_posw = 1.0, prob_dposom = 1.0)
single_gene_postpred(emp_folder, postlist_name, "uninfm2a", exp_dir + "simu30_shrink_posw01", prop_pos = 0.30, burnin=100, shrink_posw = 0.1, shrink_dposom = 1.0, prob_posw = 1.0, prob_dposom = 1.0)

# multi gene post pred: post pred itself already done (scripts written at step 1, should be launched after multigene runs are complete)
# this only formats the output
# simulating under hard or soft shrinkage

multi_gene_postpred(emp_folder, "sharedmm2a", exp_dir + "simushared")
multi_gene_postpred(emp_folder, "shrunkenmm2a", exp_dir + "simushrunken")

# prepare subexperiment folders

for simu in ["simu30", "simu10", "simu03"]:
    # prepare_subexperiment(exp_dir + simu, prep_single = True)
    prepare_subexperiment(exp_dir + simu, prep_single = False)

for simu in ["simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01"]:
    prepare_subexperiment(exp_dir + simu, prep_single = False)

for simu in ["simushared", "simushrunken"]
    prepare_subexperiment(exp_dir + simu, prep_single = False)

