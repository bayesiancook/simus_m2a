#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
res_dir = exp_dir + "results/"
with_sites = True

simu = "simu30"

# m2a_postanalysis(exp_dir + simu, ["uninfm2a"], ["indmm2a"], dlnlmin = 4.89, genepp_cutoff = 0.5, sitepp_cutoff = 0.9, outname = res_dir + "siteb" + simu, with_sites = with_sites)

m2a_postanalysis(exp_dir + simu, ["infm2a"], [], dlnlmin = 4.89, genepp_cutoff = 0.5, sitepp_cutoff = 0.9, outname = res_dir + "site90infm2a" + simu, with_sites = with_sites, refname = "infm2a")
#m2a_postanalysis(exp_dir + simu, [], ["indmm2a"], dlnlmin = 4.89, dlnlmax = 10.0, genepp_cutoff = 0.5, sitepp_cutoff = 0.9, outname = res_dir + "margsite90" + simu, with_sites = with_sites)


