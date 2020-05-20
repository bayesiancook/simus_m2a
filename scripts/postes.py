#! /usr/bin/python3.5

import sys
import numpy
import os

import parsecodeml
import parsem2a
import parsemm2a
import simuparams

def m2a_postes(exp_folder, single_basename, multi_basename, outname = "m2a_postanalysis", single_burnin = 100, multi_burnin = 500, esmode = "excess", cutoff_list = default_cutoff_list):

    exp_dir = exp_folder + "/"
    codeml_dir = exp_dir + "codeml/"
    single_dir = exp_dir + "singlegene/"
    multi_dir = exp_dir + "multigene/"

    listname = "all.list"

    single_result = dict()

    # get gene list
    with open(exp_dir + listname, 'r') as listfile:
        genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

    ngene = len(genelist)

    # get number of taxa and sites per gene
    gene_ntaxa = dict()
    gene_nsite = dict()
    fromsimu = True
    for gene in genelist:
        with open(exp_dir + "data/" + gene + ".ali") as genefile:
            [ntaxa, nsite] = genefile.readline().rstrip('\n').split()
            gene_ntaxa[gene] = int(ntaxa)
            gene_nsite[gene] = int(nsite)
        if not os.path.exists(exp_dir + "data/" + gene + ".trueparam"):
            fromsimu = False
        if not os.path.exists(exp_dir + "data/" + gene + ".truesiteom"):
            fromsimu = False

    # get true param values (if from simu)
    # and calculate empirical means and variances from these true values
    if fromsimu:
        print("data are from simulation, getting true param values..")
        [truepurw, trueposw, truepurom, truedposom, truesiteom, truepos2] = simuparams.get_true_params(exp_folder)
        trueposom = { gene : 1 + dposom for (gene,dposom) in truedposom.items()}

        trueexcess = {gene : trueposw[gene] * dposom for (gene,dposom) in truedposom.items()}

    score = dict()
    posw = dict()
    posom = dict()
    excess = dict()

    namelist = ["codeml"] + single_basename + multi_basename;

    print("parsing codeml")
    name = "codeml"
    codeml_res = parsecodeml.parse_list(codeml_dir, genelist)
    score[name] = codeml_res[0]
    posw[name] = codeml_res[1]
    posom[name] = codeml_res[2]
    excess[name] = codeml_res[7]

    print("parsing single gene analyses")
    for name in single_basename:
        print(name)
        single_res = parsem2a.parse_list(name, genelist, single_burnin, path=single_dir)
        score[name] = single_res[0]
        posw[name] = single_res[1]
        posom[name] = single_res[2]
        excess[name] = single_res[10]

    print("parsing multi gene analyses")
    for name in multi_basename:
        print(name)
        multi_res = parsemm2a.parse_list(name, multi_burnin, path=multi_dir, with_sites = False)
        score[name] = multi_res[0]
        posw[name] = multi_res[1]
        posom[name] = multi_res[2]
        excess[name] = multi_res[10]

    print("gene-level es")

    gene_ndisc = dict()
    gene_cumules = dict()
    gene_truecumules = dict()
    gene_trueesfrac = dict()
    gene_efdr = dict()
    gene_fdr = dict()

    for name in namelist:
        if esmode == "pos":
            [gene_ndisc[name], gene_cumules[name], gene_truecumules[name], gene_trueesfrac[name], gene_efdr[name], gene_fdr[name]] = gene_es(score[name], posw[name], trueposw, gene_nsite, outname + "_" + name, cutoff_list = cutoff_list)
        else:
            [gene_ndisc[name], gene_cumules[name], gene_truecumules[name], gene_trueesfrac[name], gene_efdr[name], gene_fdr[name]] = gene_es(score[name], excess[name], trueexcess, gene_nsite, outname + "_" + name, cutoff_list = cutoff_list)

    return [gene_ndisc, gene_cumules, gene_truecumules, gene_trueesfrac, gene_efdr, gene_fdr]


def full_m2a_postes(exp_folder, simu_list, single_basename, multi_basename, outname, single_burnin = 100, multi_burnin = 500, with_tex = False, esmode = "excess", cutoff_list = default_cutoff_list):

    exp_dir = exp_folder + "/"
    res_dir = exp_dir + "results/"

    simu_ndisc = dict()
    simu_cumules = dict()
    simu_truecumules = dict()
    simu_trueesfrac = dict()
    simu_efdr = dict()
    simu_fdr = dict()

    for simu in simu_list:
        print(simu)
        [simu_ndisc[simu], simu_cumules[simu], simu_truecumules[simu], simu_trueesfrac[simu], simu_efdr[simu], simu_fdr[simu]] = m2a_postes(exp_dir + simu, single_basename, multi_basename, single_burnin = single_burnin, multi_burnin = multi_burnin, outname = res_dir + simu, esmode = esmode, cutoff_list = cutoff_list)

    namelist = single_basename + multi_basename
    # namelist = ["codeml"] + single_basename + multi_basename

    with open(res_dir + outname + ".es_summary", 'w') as outfile:

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            outfile.write('{0:^24.2f}'.format(1-cutoff))
        outfile.write("\n")
        outfile.write("\n")

        for simu in simu_list:

            outfile.write("{0:>18s}".format(simu))
            for cutoff in cutoff_list:
                outfile.write(" {0:>5s}".format("n"))
                outfile.write(" {0:>5s}".format("%es"))
                outfile.write(" {0:>5s}".format("efdr"))
                outfile.write(" {0:>5s}".format("fdr"))
            outfile.write("\n")

            for name in namelist:

                outfile.write("{0:>18s}".format(name))

                for cutoff in cutoff_list:
                    outfile.write(" {0:5d}".format(simu_ndisc[simu][name][cutoff]))
                    outfile.write(" {0:5.2f}".format(simu_trueesfrac[simu][name][cutoff]))
                    outfile.write(" {0:5.2f}".format(simu_efdr[simu][name][cutoff]))
                    outfile.write(" {0:5.2f}".format(simu_fdr[simu][name][cutoff]))
                outfile.write("\n")
            outfile.write("\n")

