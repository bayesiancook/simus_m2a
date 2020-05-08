#! /usr/bin/python3.5

import sys
import os

from postanalysis import * 

exp_folder = sys.argv[1]
exp_dir = exp_folder + "/"
res_dir = exp_dir + "results/"
with_sites = False

fdr_cutoff_list = [0.05, 0.1, 0.3, 0.5]

# a restricted post analysis under codeml and multi gene models

multi_basename = ["unconsindmm2a", "indmm2a", "unconsshrunkenmm2a", "shrunkenmm2a", "sharedmm2a"]
namelist = ["df1_codeml", "df2_codeml", "mixdf1_codeml"] + multi_basename

simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_dposom03", "simu30_shrink_dposom01", "simu30_shrink_posw03", "simu30_shrink_posw01", "simushrunken", "simushared"]

simu_ndisc = dict()
simu_fdr = dict()

for simu in simu_list:
    print(simu)
    [simu_ndisc[simu], simu_fdr[simu]] = m2a_postanalysis(exp_dir + simu, [], multi_basename, outname = res_dir + "multi" + simu, with_sites = with_sites, cutoff_list = fdr_cutoff_list)


with open(res_dir + "fdrtable.tex", 'w') as outfile:

    tabul = 'rr' + 'c' * 2 * len(fdr_cutoff_list)
    outfile.write("\\begin{{tabular}}{{{0}}}\n".format(tabul))

    outfile.write("& & \\multicolumn{{ {0} }}{{c}}{{target FDR}}".format(2*len(fdr_cutoff_list)))
    outfile.write(r'\\')
    outfile.write("\n")

    outfile.write(r'&')
    for cutoff in fdr_cutoff_list:
        outfile.write("& \\multicolumn{{2}}{{c}}{{${0:5.2f}$}}".format(cutoff))
    outfile.write(r'\\')
    outfile.write("\n")

    outfile.write(r'simulation & method ')
    for cutoff in fdr_cutoff_list:
        outfile.write(r'& nd & fdr ')
    outfile.write(r'\\')
    outfile.write("\n")

    outfile.write("\\hline\n")

    for simu in simu_list:

        for name in namelist:

            if name[-6:] == "codeml":
                outfile.write("{0:>18s}".format(simu.replace("_", "\\_")))

            outfile.write("&{0:>18s}".format(name))
            for cutoff in fdr_cutoff_list:
                ndisc = simu_ndisc[simu][name][cutoff]
                if ndisc:
                    fdr = simu_fdr[simu][name][cutoff]
                    outfile.write("& ${0:5d}$ & ${1:5.2f}$ ".format(ndisc, fdr))

            outfile.write(r'\\')
            outfile.write("\n")

        outfile.write("\\hline\n")

    outfile.write(r'\end{tabular}{}')
    outfile.write("\n")


