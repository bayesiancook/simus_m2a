#! /usr/bin/python3.5

import postpredsimu
import sys

exp_folder = sys.argv[1]
listname = sys.argv[2]
lnl_min = float(sys.argv[3])
shrink_posw = float(sys.argv[4])
shrink_dposom = float(sys.argv[5])
# lnl_max = float(sys.argv[4])
target_name = sys.argv[6]

basename = "uninfm2a"

postpredsimu.single_gene_postpred(exp_folder, listname, basename, target_name, lnl_min = lnl_min, shrink_posw = shrink_posw, shrink_dposom = shrink_dposom, burnin=50)
# postpredsimu.single_gene_postpred(exp_folder, listname, basename, target_name, lnl_min = lnl_min, lnl_max = lnl_max, burnin=50)

