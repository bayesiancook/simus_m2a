#! /usr/bin/python3.5

import sys
import numpy
import os

import parsecodeml
import parsem2a
import parsemm2a
import simuparams

cutoff_list = [0.05, 0.1, 0.3, 0.5]

def gene_oma(score, oma, posw, posom, trueoma, trueposw, trueposom, nsite, outname):

    gene_ndisc = dict()
    gene_truefrac = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_truefrac[cutoff] = 0

    ngene = len(oma)
    fromsimu = len(trueoma)

    with open(outname + ".sortedoma", 'w') as outfile:

        if fromsimu:

            outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}  {7:5s}  {8:5s}  {9:5s}  {10:5s}  {11:5s}\n".format("#gene", "nsite", "e-oma", "oma", "e-%oma", "%oma", "%genes", "pp", "e-p+", "p+", "-om+", "om+"))
            totnsite = sum([ns for (gene,ns) in nsite.items()])
            grandtotoma = sum([oma*nsite[gene] for (gene,oma) in oma.items()])
            grandtottrueoma = sum([oma*nsite[gene] for (gene,oma) in trueoma.items()])
            tottruea = 0
            tota = 0
            n = 0

            for (gene,oma) in sorted(oma.items(), key=lambda kv: kv[1], reverse=True):

                a = oma * nsite[gene]
                tota = tota + a
                truea = trueoma[gene] * nsite[gene]
                tottruea = tottruea + truea
                n = n + 1
                fracoma = tota/grandtotoma
                fractrueoma = tottruea/grandtottrueoma

                outfile.write("{0:18s}  {1:5d}  {2:5.3f}  {3:5.3f}  {4:5.2f}  {5:5.2f}  {6:5.2f}  {7:5.3f}  {8:5.3f}  {9:5.3f}  {10:5.3f}  {11:5.3f}\n".format(
                    gene, nsite[gene], oma, trueoma[gene], fracoma, fractrueoma, n/ngene, score[gene], posw[gene], trueposw[gene], posom[gene], trueposom[gene]))

                for cutoff in cutoff_list:
                    if fracoma < 1-cutoff:
                        gene_ndisc[cutoff] = n / ngene
                        gene_truefrac[cutoff] = fractrueoma

        # else:
            # outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}  {7:5s}\n".format("gene", "nsite", "om_a", "%oma", "%genes", "pp", "p+", "om+"))

    return [gene_ndisc, gene_truefrac]




def m2a_postoma(exp_folder, single_basename, multi_basename, outname = "m2a_postanalysis", single_burnin = 100, multi_burnin = 500):

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
        trueoma = {gene : trueposw[gene] * dposom for (gene,dposom) in truedposom.items()}


    score = dict()
    posw = dict()
    posom = dict()
    oma = dict()

    namelist = ["codeml"] + single_basename + multi_basename;

    print("parsing codeml")
    name = "codeml"
    codeml_res = parsecodeml.parse_list(codeml_dir, genelist)
    score[name] = codeml_res[0]
    posw[name] = codeml_res[1]
    posom[name] = codeml_res[2]
    oma[name] = codeml_res[7]

    print("parsing single gene analyses")
    for name in single_basename:
        print(name)
        single_res = parsem2a.parse_list(name, genelist, single_burnin, path=single_dir)
        score[name] = single_res[0]
        posw[name] = single_res[1]
        posom[name] = single_res[2]
        oma[name] = single_res[10]

    print("parsing multi gene analyses")
    for name in multi_basename:
        print(name)
        multi_res = parsemm2a.parse_list(name, multi_burnin, path=multi_dir, with_sites = False)
        score[name] = multi_res[0]
        posw[name] = multi_res[1]
        posom[name] = multi_res[2]
        oma[name] = multi_res[10]

    print("gene-level oma")

    gene_ndisc = dict()
    gene_truefrac = dict()

    for name in namelist:
        [gene_ndisc[name], gene_truefrac[name]] = gene_oma(score[name], oma[name], posw[name], posom[name], trueoma, trueposw, trueposom, gene_nsite, outname + "_" + name)

    return [gene_ndisc, gene_truefrac]


def full_m2a_postoma(exp_folder, simu_list, single_basename, multi_basename, outname, single_burnin = 100, multi_burnin = 500, with_tex = False):

    exp_dir = exp_folder + "/"
    res_dir = exp_dir + "results/"

    simu_ndisc = dict()
    simu_truefrac = dict()

    for simu in simu_list:
        print(simu)
        [simu_ndisc[simu], simu_truefrac[simu]] = m2a_postoma(exp_dir + simu, single_basename, multi_basename, single_burnin = single_burnin, multi_burnin = multi_burnin, outname = res_dir + simu)

    namelist = ["codeml"] + single_basename + multi_basename

    with open(res_dir + outname + ".oma_summary", 'w') as outfile:

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            outfile.write('{0:^12.2f}'.format(1-cutoff))
        outfile.write("\n")
        outfile.write("\n")

        for simu in simu_list:

            outfile.write("{0:>18s}".format(simu))
            for cutoff in cutoff_list:
                outfile.write(" {0:>5s}".format("%gene"))
                outfile.write(" {0:>5s}".format("%oma"))
            outfile.write("\n")

            for name in namelist:

                outfile.write("{0:>18s}".format(name))

                for cutoff in cutoff_list:
                    outfile.write(" {0:5.2f}".format(simu_ndisc[simu][name][cutoff]))
                    outfile.write(" {0:5.2f}".format(simu_truefrac[simu][name][cutoff]))
                outfile.write("\n")
            outfile.write("\n")

