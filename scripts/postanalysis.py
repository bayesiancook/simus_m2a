#! /usr/bin/python3.5

import sys
import numpy
import os
import parsecodeml
import parsem2a
import parsemm2a
import simuparams
from scipy.stats import chi2

def gene_codeml_fdr(score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    with open(outname + "_codeml.genefdr", 'w') as outfile:

        # header 
        if fromsimu:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\t{6:7s}\n".format("#gene", "dlnL", "pval", "e-fdr", "fdr", "tpr", "fpr"))
        else:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\n".format("gene", "dlnL", "pval", "e-fdr"))

        n = 0
        fp = 0
        tp = 0
        for (gene,dlnl) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):

            if dlnl>0:
                n = n + 1
                pval = (1 - chi2.cdf(2*dlnl,1))/2
                efdr = ngene*pval/n
                if efdr > 1:
                    efdr = 1
                if fromsimu:
                    if truepos[gene] > 0:
                        tp = tp + 1
                    else:
                        fp = fp + 1
                    # false discovery rate
                    fdr = fp / n
                    # sensitivity or recall
                    tpr = tp / ntrue
                    # 1 - specificity: false positive rate
                    fpr = fp / nfalse
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\t{6:7.5f}\n".format(gene,dlnl,pval,efdr,fdr,tpr,fpr))
                else:
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,dlnl,pval,efdr))

def gene_fdr(cutoff_list, score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    gene_ndisc = dict()
    gene_fp = dict()
    gene_efdr = dict()
    gene_etpr = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fp[cutoff] = 0
        gene_efdr[cutoff] = 0
        gene_etpr[cutoff] = 0

    with open(outname + ".genefdr", 'w') as outfile:

        # header
        if fromsimu:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\t{6:7s}\n".format("#gene", "pp", "e-fdr", "e-tpr", "fdr", "tpr", "fpr"))
        else:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\t{3:7s}\n".format("gene", "pp", "e-fdr", "e-tpr"))

        # estimated total number of positively selected genes
        etotp = sum([pp for (gene,pp) in score.items()])

        fp = 0
        tp = 0
        n = 0
        totpp = 0

        for (gene,pp) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):
            n = n + 1
            # total post prob for positive selection thus far: estimated number of true positives in selected set
            totpp = totpp + pp
            # estimated rate of false discovery: 1 - FP / N
            efdr = 1 - totpp/n
            # estimated sensitivity (true positive rate): TP / TOT NUMBER OF POS SEL GENES
            etpr = totpp / etotp

            if fromsimu:
                if truepos[gene] > 0:
                    tp = tp + 1
                else:
                    fp = fp + 1
                # false discovery rate
                fdr = fp / n
                # sensitivity or recall
                tpr = tp / ntrue
                # 1 - specificity: false positive rate
                fpr = fp / (ngene - ntrue)
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\t{6:7.5f}\n".format(gene,score[gene],efdr,etpr,fdr,tpr,fpr))
            else:
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,score[gene],efdr,etpr))

            for cutoff in cutoff_list:
                if pp > cutoff:
                    gene_ndisc[cutoff] = n
                    gene_fp[cutoff] = fp
                    gene_efdr[cutoff] = efdr
                    gene_etpr[cutoff] = etpr

    return [gene_ndisc, gene_fp, gene_efdr, gene_etpr]


def method_gene_fdr(cutoff_list, namelist, score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    gene_ndisc = dict()
    gene_fp = dict()
    gene_efdr = dict()
    gene_etpr = dict()

    for name in namelist:
        [gene_ndisc[name], gene_fp[name], gene_efdr[name], gene_etpr[name]] = gene_fdr(cutoff_list, score[name], truepos, outname + "_" + name);
    with open(outname + ".genefdr", 'w') as outfile:

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            if fromsimu:
                outfile.write("{0:^29.2f}".format(cutoff))
            else:
                outfile.write("{0:^19.2f}".format(cutoff))
        outfile.write("\n")

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            if fromsimu:
                outfile.write("  {0:>5s} {1:>5s} {2:>5s} {3:>5s} {4:>5s}".format("disc", "efdr", "fdr", "etpr", "tpr"))
            else:
                outfile.write("  {0:>5s} {1:>5s} {2:>5s}".format("disc", "efdr", "etpr"))
        outfile.write("\n")

        for name in namelist:

            if name != "codeml":
                outfile.write("{0:>18s}".format(name))

                for cutoff in cutoff_list:
                    ndisc = gene_ndisc[name][cutoff]
                    if ndisc:
                        efdr = gene_efdr[name][cutoff]
                        etpr = gene_etpr[name][cutoff]
                        if fromsimu:
                            fdr = gene_fp[name][cutoff] / ndisc
                            tpr = (ndisc - gene_fp[name][cutoff]) / ntrue
                            outfile.write("  {0:5d} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f}".format(ndisc, efdr, fdr, etpr, tpr))
                        else:
                            outfile.write("  {0:5d} {1:5.2f} {2:5.2f}".format(ndisc, efdr, etpr))

            outfile.write("\n")

def m2a_postanalysis(exp_folder, single_basename, multi_basename, outname = "m2a_postanalysis", burnin = 20, dlnlmin = 0, min_omega = 1.0):

    cutoff_list = [0.1,0.3, 0.5, 0.7, 0.9]

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
        (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = simuparams.get_empirical_hyperparams(truepurw, trueposw, truepurom, truedposom)

    # for each method (codeml, m2a, mm2a, under any prior or settings)
    # and for all genes
    # we want posw and pp, posom and min max, and selected sites
    # posw[method_name][gene_name] ...

    posw = dict()
    score = dict()
    score2 = dict()
    score3 = dict()
    posom = dict()
    minposom = dict()
    maxposom = dict()
    selectedsites = dict()
    sitepp = dict()
    hyperparams = dict()

    namelist = []

    print("parsing codeml")
    # parsing codeml results
    name = "codeml"
    codeml_res = parsecodeml.parse_list(codeml_dir, genelist)
    [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name]] = codeml_res[0:7]
    namelist.append(name)

    print("parsing single gene analyses")
    # parsing m2a results
    for name in single_basename:
        print(name)
        single_res = parsem2a.parse_list(name, genelist, int(burnin), path=single_dir, min_omega = min_omega)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name]] = single_res[0:9]
        namelist.append(name)

    print("parsing multi gene analyses")
    # parsing mm2a results
    for name in multi_basename:
        print(name)
        multi_res = parsemm2a.parse_list(name, int(burnin), path=multi_dir)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name], hyperparams[name]] = multi_res[0:10]
        namelist.append(name)

    print("printing out sorted list in .sortedparams")

    codeml_dlnl = score["codeml"]

    with open(outname + ".sortedparams", 'w') as geneoutfile:

        # header 

        geneoutfile.write("methods:\n")
        for i,name in enumerate(namelist):
            geneoutfile.write("{}\t{}\n".format(i,name))
        geneoutfile.write("\n")

        geneoutfile.write("{:15}".format("gene"))
        geneoutfile.write("  {:>6}".format("dlnl"))
        for i in range(len(namelist)-1):
            geneoutfile.write("{:>7}{:1}".format("pp_",i+1))

        geneoutfile.write(" ")

        if fromsimu:
            geneoutfile.write("  {:>5}".format("w*"))

        for i in range(len(namelist)):
            geneoutfile.write("  {:>4}{:1}".format("w_",i))

        geneoutfile.write(" ")

        if fromsimu:
            geneoutfile.write("  {:>6}".format("om*"))

        for i in range(len(namelist)):
            geneoutfile.write("  {:>5}{:1}".format("om_",i))

        geneoutfile.write(" ")

        geneoutfile.write("\n")
        geneoutfile.write("\n")

        #
        # gene list in decreasing order of dlnl
        #

        for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):

            geneoutfile.write("{0:15s}".format(gene))
            for name in namelist:
                geneoutfile.write("  {0:6.2f}".format(score[name][gene]))

            geneoutfile.write(" ")

            if fromsimu:
                geneoutfile.write("  {0:5.3f}".format(trueposw[gene]))

            for name in namelist:
                geneoutfile.write("  {0:5.3f}".format(posw[name][gene]))

            geneoutfile.write(" ")

            if fromsimu:
                geneoutfile.write("  {0:6.2f}".format(1.0 + truedposom[gene]))

            for name in namelist:
                geneoutfile.write("  {0:6.2f}".format(posom[name][gene]))

            # site-level statistics -- deactivated

            #for cutoff in cutoff_list:
            #    geneoutfile.write(" ")
            #    for name in namelist:
            #        geneoutfile.write("  {0:2d}".format(nsel[cutoff][name][gene]))
            #        if fromsimu:
            #            if name == "codeml":
            #                geneoutfile.write(" ({0:2d})".format(tdr[cutoff][name][gene]))
            #            else:
            #                geneoutfile.write(" ({0:2d} ; {1:3.1f})".format(tdr[cutoff][name][gene], etdr[cutoff][name][gene]))
            geneoutfile.write("\n")

    #with open(outname + ".sitepp", 'w') as siteoutfile:

        #
        # gene list in decreasing order of dlnl
        #

    #    for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):
    #        if lnl >= dlnlmin:
    #            for (i,pp) in selectedsites["codeml"][gene].items():
    #                siteoutfile.write("{0}\t{1}\t{2}".format(gene,i+1,pp))
    #                for name in namelist:
    #                    if name != "codeml":
    #                        siteoutfile.write("\t{0}".format(sitepp[name][gene][i]))
    #                siteoutfile.write("\n")

    hypernamelist = ['purom_mean', 'purom_invconc', 'dposom_mean', 'dposom_invshape', 'purw_mean', 'purw_invconc', 'posw_mean', 'posw_invconc', 'pi']

    print("printing out estimated hyperparameters")
    with open(outname + ".hyper", 'w') as outfile:

        # posterior mean and credible interval for hyperparams of multigene runs
        outfile.write("{0:20s}".format("param"))
        for name in hyperparams:
            outfile.write("\t{0:^24s}".format(name))
        outfile.write('\n')
        for i,name in enumerate(hypernamelist):
            outfile.write("{0:20s}".format(name))
            for name in hyperparams:
                (mean,min,max) = hyperparams[name][i]
                outfile.write("\t{0:6.4f} ({1:6.4f} , {2:6.4f})".format(mean,min,max))
            outfile.write('\n')

        outfile.write('\n')

        # mean and stdev of distribution of true versus post mean estimates for the 4 key parameters
        # outfile.write("{0:10s}\t{1:^6s}\t{2:^6s}\t{3:^6s}\t{4:^6s}\t{5:^6s}\n".format("method", "pi", "posw", "stdev", "posom", "stedv"))
        # if fromsimu:
        #    outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format("simu", pi, posw_mean, numpy.sqrt(posw_var), 1.0 + dposom_mean, numpy.sqrt(dposom_var)))

        #for name in namelist:
        #    (est_pi, est_purw_mean, est_purw_var, est_purw_invconc, est_posw_mean, est_posw_var, est_posw_invconc, est_purom_mean, est_purom_var, est_purom_invconc, est_posom_mean, est_posom_var, est_posom_invshape) = simuparams.get_empirical_hyperparams(posw[name], posw[name], posom[name], posom[name])
        #    outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format(name, est_pi, est_posw_mean, numpy.sqrt(est_posw_var), est_posom_mean, numpy.sqrt(est_posom_var)))

    print("gene-level fdr")
    truepos = []
    if fromsimu:
        truepos = trueposw

    gene_codeml_fdr(score["codeml"], truepos, outname)
    namelist2 = [name for name in namelist if name != "codeml"]
    method_gene_fdr(cutoff_list, namelist2, score, truepos, outname)

    print("done")

