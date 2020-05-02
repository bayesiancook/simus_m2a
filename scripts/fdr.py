#! /usr/bin/python3.5

import sys
import numpy
import os
from scipy.stats import chi2

def gene_codeml_fdr(cutoff_list, score, truepos, outname):

    gene_ndisc = dict()
    gene_fp = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fp[cutoff] = 0

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

                    for cutoff in cutoff_list:
                        if efdr < cutoff:
                            gene_ndisc[cutoff] = n
                            gene_fp[cutoff] = fp

                else:
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,dlnl,pval,efdr))

    return [gene_ndisc, gene_fp]

def gene_fdr(cutoff_list, score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    gene_ndisc = dict()
    gene_fp = dict()
    gene_etpr = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fp[cutoff] = 0
        gene_etpr[cutoff] = 0

    with open(outname + ".genefdr", 'w') as outfile:

        # header
        if fromsimu:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\n".format("#gene", "pp", "e-tpr", "fdr", "tpr", "fpr"))
        else:
            outfile.write("{0:15s}\t{1:7s}\t{2:7s}\n".format("gene", "pp", "e-tpr"))

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
                if efdr < cutoff:
                    gene_ndisc[cutoff] = n
                    gene_fp[cutoff] = fp
                    gene_etpr[cutoff] = etpr

    return [gene_ndisc, gene_fp, gene_etpr]


def method_gene_fdr(cutoff_list, namelist, score, truepos, outname):

    ngene = len(score)
    fromsimu = len(truepos)
    ntrue = len([x for (gene,x) in truepos.items() if x > 0])
    nfalse = ngene - ntrue

    gene_ndisc = dict()
    gene_fp = dict()
    gene_etpr = dict()

    for name in namelist:
        if name == "codeml":
            [gene_ndisc[name], gene_fp[name]] = gene_codeml_fdr(cutoff_list, score[name], truepos, outname + "_" + name);
        else:
            [gene_ndisc[name], gene_fp[name], gene_etpr[name]] = gene_fdr(cutoff_list, score[name], truepos, outname + "_" + name);

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
                outfile.write("  {0:>5s} {1:>5s} {2:>5s} {3:>5s}".format("disc", "fdr", "etpr", "tpr"))
            else:
                outfile.write("  {0:>5s} {1:>5s}".format("disc", "etpr"))
        outfile.write("\n")

        for name in namelist:

            outfile.write("{0:>18s}".format(name))

            for cutoff in cutoff_list:
                ndisc = gene_ndisc[name][cutoff]
                if ndisc:
                    if fromsimu:
                        fdr = gene_fp[name][cutoff] / ndisc
                        tpr = (ndisc - gene_fp[name][cutoff]) / ntrue
                        if name == "codeml":
                            outfile.write("  {0:5d} {1:5.2f} {2:^5s} {3:5.2f}".format(ndisc, fdr, "-", tpr))
                        else:
                            etpr = gene_etpr[name][cutoff]
                            outfile.write("  {0:5d} {1:5.2f} {2:5.2f} {3:5.2f}".format(ndisc, fdr, etpr, tpr))

                    else:
                        if name == "codeml":
                            outfile.write("  {0:5d} {1:^5s}".format(ndisc, "-"))
                        else:
                            etpr = gene_etpr[name][cutoff]
                            outfile.write("  {0:5d} {1:5.2f}".format(ndisc, etpr))

            outfile.write("\n")
