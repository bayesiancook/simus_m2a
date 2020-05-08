#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
res_dir = exp_dir + "results/"
with_sites = True

simu = "simu30"

#m2a_postanalysis(exp_dir + simu, [], ["indmm2a"], dlnlmin = 4.89, genepp_cutoff = 0.5, sitepp_cutoff = 0.9, outname = res_dir + "site90" + simu, with_sites = with_sites)
#m2a_postanalysis(exp_dir + simu, [], ["indmm2a"], dlnlmin = 4.89, genepp_cutoff = 0.5, sitepp_cutoff = 0.95, outname = res_dir + "site95" + simu, with_sites = with_sites)

m2a_postanalysis(exp_dir + simu, [], ["indmm2a"], dlnlmin = 4.89, dlnlmax = 10.0, genepp_cutoff = 0.5, sitepp_cutoff = 0.9, outname = res_dir + "margsite90" + simu, with_sites = with_sites)
m2a_postanalysis(exp_dir + simu, [], ["indmm2a"], dlnlmin = 4.89, dlnlmax = 10.0, genepp_cutoff = 0.5, sitepp_cutoff = 0.95, outname = res_dir + "margsite95" + simu, with_sites = with_sites)

#m2a_postanalysis(exp_dir + simu, [], ["unconsindmm2a"], dlnlmin = 4.89, genepp_cutoff = 0.589, refname = "unconsindmm2a", sitepp_cutoff = 0.9, outname = res_dir + "unconssite90" + simu, with_sites = with_sites)
#m2a_postanalysis(exp_dir + simu, [], ["unconsindmm2a"], dlnlmin = 4.89, genepp_cutoff = 0.589, refname = "unconsindmm2a", sitepp_cutoff = 0.95, outname = res_dir + "unconssite95" + simu, with_sites = with_sites)
