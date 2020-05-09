#! /usr/bin/python3.5

#
# post analysis of experiments on original empirical data
# 

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
emp_folder = exp_folder + "/empirical"
res_folder = exp_folder + "/results"

fdr_cutoff_list = [0.05, 0.1, 0.3, 0.5]

if not os.path.exists(res_folder):
    print("creating results folder")
    os.system("mkdir " + res_folder)

# post analysis

single_basename = ["uninfm2a", "uninfpi50m2a", "uninfpi10m2a", "uninfpi02m2a", "subjpi50m2a", "subjpi10m2a", "subjpi02m2a"]
multi_basename = ["indmm2a", "shrunkenmm2a", "sharedmm2a", "unconsindmm2a", "unconsshrunkenmm2a"]

[method_ndisc, method_fdr, method_efnr, method_fnr] = m2a_postanalysis(emp_folder, single_basename, multi_basename, outname = exp_folder + "/results/emp")

namelist = ["df1_codeml", "df2_codeml", "mixdf1_codeml"] + single_basename + multi_basename

with open(res_folder + "/emp.summary", 'w') as outfile:

    outfile.write("{0:>18s}".format(""))
    for cutoff in fdr_cutoff_list:
        outfile.write("{0:^14.2f}".format(cutoff))
    outfile.write("\n")

    outfile.write("{0:>18s}".format(""))
    for cutoff in fdr_cutoff_list:
        outfile.write("  {0:>5s}  {1:>5s}".format("disc", "e-fnr"))
    outfile.write("\n")

    for name in namelist:

        outfile.write("{0:>18s}".format(name))

        for cutoff in fdr_cutoff_list:
            if name[-6:] == "codeml":
                outfile.write("  {0:5d}  {1:^5s}".format(method_ndisc[name][cutoff], "-"))
            else:
                outfile.write("  {0:5d}  {1:5.2f}".format(method_ndisc[name][cutoff], 
                    method_efnr[name][cutoff]))

        outfile.write("\n")

