#! /usr/bin/python3.5

#
# post pred on multi gene runs: scripts already produced at step1
# here: assume that those scripts have already been run but still has to format the output
# 
# then: prepare the subexperiment folders 
# 

import sys
import os

from postpredsimu import *
from prepare_subexperiment import *

exp_folder = sys.argv[1]

exp_dir = exp_folder + "/"

# multi gene post pred: post pred itself already done (scripts written at step 1, should be launched after multigene runs are complete)
# this only formats the output
# simulating under hard or soft shrinkage

multi_gene_postpred(emp_folder, "indmm2a", exp_dir + "simuind")
multi_gene_postpred(emp_folder, "sharedmm2a", exp_dir + "simushared")
multi_gene_postpred(emp_folder, "shrunkenmm2a", exp_dir + "simushrunken")

for simu in ["simushared", "simushrunken", "simuind"]:
    prepare_subexperiment(exp_dir + simu, prep_single = False)

