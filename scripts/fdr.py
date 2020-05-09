#! /usr/bin/python3.5

import sys
import numpy
import os
from scipy.stats import chi2

# based on list of dlnl from codeml, compute p-value and apply BH algorithm
# returns, for a series of target FDR thresholds:
# - number of discoveries, 
# - false discovery rate
# - false negative rate
# 
# assumed null distribution for computing p-values:
# - chi_mode = df1    : chi2 with 1 df
# - chi_mode = df2    : chi2 with 2 df
# - chi_mode = mixdf1 : mix of 0.5 point mass at 0 and 0.5 chi2 with 1df

def gene_codeml_fdr(cutoff_list, score, truepos, outname, chi_mode = "df2"):

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_fnr = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fdr[cutoff] = 0
        gene_fnr[cutoff] = 0

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    with open(outname + "_codeml.genefdr", 'w') as outfile:

        # header 
        if fromsimu:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\n".format("#gene", "dlnL", "pval", "e-fdr", "fdr", "fnr"))
        else:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\n".format("gene", "dlnL", "pval", "e-fdr"))

        n = 0
        fp = 0
        tp = 0
        fdr = 0

        for (gene,dlnl) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):

            if dlnl>0:
                n = n + 1
                pval = 0
                if chi_mode == "df2":
                    pval = 1 - chi2.cdf(2*dlnl,2)
                else:
                    if chi_mode == "df1":
                        pval = 1 - chi2.cdf(2*dlnl,1)
                    else:
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
                    # false negative rate
                    fnr = 1 - tpr
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\n".format(gene,dlnl,pval,efdr,fdr,fnr))

                    for cutoff in cutoff_list:
                        if efdr < cutoff:
                            gene_ndisc[cutoff] = n
                            gene_fdr[cutoff] = fdr
                            gene_fnr[cutoff] = fnr

                else:
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,dlnl,pval,efdr))

                    for cutoff in cutoff_list:
                        if efdr < cutoff:
                            gene_ndisc[cutoff] = n

    return [gene_ndisc, gene_fdr, gene_fnr]

# based on list of pp from a bayesian analysis,
# compute bayes e-fdr (expected or advertised FDR)
# returns, for a series of target FDR thresholds:
# - number of discoveries, 
# - fdr
# - estimated fnr
# - fnr

def gene_bayes_fdr(cutoff_list, score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_efnr = dict()
    gene_fnr = dict()

    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fdr[cutoff] = 0
        gene_efnr[cutoff] = 0
        gene_fnr[cutoff] = 0

    with open(outname + ".genefdr", 'w') as outfile:

        # header
        if fromsimu:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\n".format("#gene", "pp", "e-fdr", "fdr", "e-fnr", "fnr"))
        else:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\t{3:7s}\n".format("gene", "pp", "e-fdr", "e-fnr"))

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
            # estimated false negative rate
            efnr = 1 - etpr

            fdr = 0
            fnr = 0

            if fromsimu:
                if truepos[gene] > 0:
                    tp = tp + 1
                else:
                    fp = fp + 1
                # false discovery rate
                fdr = fp / n
                # sensitivity or recall
                tpr = tp / ntrue
                # false negative rate
                fnr = 1 - tpr
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\n".format(gene,score[gene],efdr,fdr,efnr,fnr))
            else:
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,score[gene],efdr,efnr))

            for cutoff in cutoff_list:
                if efdr < cutoff:
                    gene_ndisc[cutoff] = n
                    gene_fdr[cutoff] = fdr
                    gene_efnr[cutoff] = efnr
                    gene_fnr[cutoff] = fnr

    return [gene_ndisc, gene_fdr, gene_efnr, gene_fnr]


# tabulate results across all methods for a given dataset

def method_gene_fdr(cutoff_list, namelist, score, truepos, outname):

    fromsimu = len(truepos)

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_efnr = dict()
    gene_fnr = dict()

    for name in namelist:
        if name == "codeml":
            for mode in ["df2", "df1", "mixdf1"]:
                [gene_ndisc[mode + "_codeml"], gene_fdr[mode + "_codeml"], gene_fnr[mode + "_codeml"]] = gene_codeml_fdr(cutoff_list, score["codeml"], truepos, outname + "_" + mode, chi_mode = mode);

        else:
            [gene_ndisc[name], gene_fdr[name], gene_efnr[name], gene_fnr[name]] = gene_bayes_fdr(cutoff_list, score[name], truepos, outname + "_" + name);

    with open(outname + ".genefdr", 'w') as outfile:

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            if fromsimu:
                outfile.write("{0:^24.2f}".format(cutoff))
            else:
                outfile.write("{0:^14.2f}".format(cutoff))
        outfile.write("\n")

        outfile.write("{0:>18s}".format(""))
        for cutoff in cutoff_list:
            if fromsimu:
                outfile.write("  {0:>5s} {1:>5s} {2:>5s} {3:>5s}".format("disc", "fdr", "efnr", "fnr"))
            else:
                outfile.write("  {0:>5s} {1:>5s}".format("disc", "efnr"))
        outfile.write("\n")

        namelist2 = ["df1_codeml", "df2_codeml", "mixdf1_codeml"] + [name for name in namelist if name != "codeml"]
        for name in namelist2:

            outfile.write("{0:>18s}".format(name))

            for cutoff in cutoff_list:
                if fromsimu:
                    if name[-6:] == "codeml":
                        outfile.write("  {0:5d} {1:5.2f} {2:^5s} {3:5.2f}".format(
                            gene_ndisc[name][cutoff], 
                            gene_fdr[name][cutoff], 
                            "-", 
                            gene_fnr[name][cutoff]))
                    else:
                        outfile.write("  {0:5d} {1:5.2f} {2:5.2f} {3:5.2f}".format(
                            gene_ndisc[name][cutoff], 
                            gene_fdr[name][cutoff], 
                            gene_efnr[name][cutoff],
                            gene_fnr[name][cutoff]))

                else:
                    if name[-6:] == "codeml":
                        outfile.write("  {0:5d} {1:^5s}".format(
                            gene_ndisc[name][cutoff],
                            "-"))
                    else:
                        outfile.write("  {0:5d} {1:5.2f}".format(
                            gene_ndisc[name][cutoff],
                            gene_efnr[name][cutoff]))

            outfile.write("\n")

    return [gene_ndisc, gene_fdr, gene_efnr, gene_fnr]

