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

def gene_codeml_fdr(cutoff_list, score, truees, gene_nsite, outname, chi_mode = "df2"):

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_fnr = dict()
    gene_truecumules = dict()
    gene_trueesfrac = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fdr[cutoff] = 0
        gene_fnr[cutoff] = 0
        gene_truecumules[cutoff] = 0
        gene_trueesfrac[cutoff] = 0

    ngene = len(score)
    fromsimu = len(truees)
    ntrue = len([x for (gene,x) in truees.items() if x > 0])
    nfalse = ngene - ntrue

    totnsite = sum([nsite for (gene,nsite) in gene_nsite.items()])
    truetotes = sum([nsite*truees[gene] for (gene,nsite) in gene_nsite.items()])

    with open(outname + "_codeml.genefdr", 'w') as outfile:

        # header 
        if fromsimu:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\t{4:7s}\t{5:7s}\n".format("#gene", "dlnL", "pval", "e-fdr", "fdr", "fnr","%es"))
        else:
            outfile.write("{0:15s}\t{1:5s}\t{2:7s}\t{3:7s}\n".format("gene", "dlnL", "pval", "e-fdr"))

        n = 0
        fp = 0
        tp = 0
        fdr = 0
        truecumules = 0
        ns = 0
        

        for (gene,dlnl) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):

            if dlnl>0:
                n = n + 1
                ns = ns + gene_nsite[gene]
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
                    if truees[gene] > 0:
                        tp = tp + 1
                    else:
                        fp = fp + 1
                    # truefrac
                    truecumules = truecumules + gene_nsite[gene]*truees[gene]
                    trueesfrac = truecumules / truetotes

                    # false discovery rate
                    fdr = fp / n
                    # sensitivity or recall
                    tpr = tp / ntrue
                    # false negative rate
                    fnr = 1 - tpr

                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\t{6:7.5f}\n".format(gene,dlnl,pval,efdr,fdr,fnr,trueesfrac))

                    for cutoff in cutoff_list:
                        if efdr < cutoff:
                            gene_ndisc[cutoff] = n
                            gene_fdr[cutoff] = fdr
                            gene_fnr[cutoff] = fnr
                            # gene_esfrac[cutoff] = esfrac
                            gene_truecumules[cutoff] = truecumules / ns * 100
                            gene_trueesfrac[cutoff] = trueesfrac

                else:
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\n".format(gene,dlnl,pval,efdr))

                    for cutoff in cutoff_list:
                        if efdr < cutoff:
                            gene_ndisc[cutoff] = n

    return [gene_ndisc, gene_fdr, gene_fnr, gene_truecumules, gene_trueesfrac]

# based on list of pp from a bayesian analysis,
# compute bayes e-fdr (expected or advertised FDR)
# returns, for a series of target FDR thresholds:
# - number of discoveries, 
# - fdr
# - estimated fnr
# - fnr

def gene_bayes_fdr(cutoff_list, score, meanes, truees, gene_nsite, outname):

    ngene = len(score)
    fromsimu = len(truees)
    ntrue = len([x for (gene,x) in truees.items() if x > 0])
    nfalse = ngene - ntrue

    totnsite = sum([nsite for (gene,nsite) in gene_nsite.items()])
    totes = sum([nsite*meanes[gene] for (gene,nsite) in gene_nsite.items()])
    truetotes = sum([nsite*truees[gene] for (gene,nsite) in gene_nsite.items()])

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_efnr = dict()
    gene_fnr = dict()
    gene_cumules = dict()
    gene_esfrac = dict()
    gene_truecumules = dict()
    gene_trueesfrac = dict()

    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_fdr[cutoff] = 0
        gene_efnr[cutoff] = 0
        gene_fnr[cutoff] = 0
        gene_cumules[cutoff] = 0
        gene_esfrac[cutoff] = 0
        gene_truecumules[cutoff] = 0
        gene_trueesfrac[cutoff] = 0

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
        cumules = 0
        truecumules = 0
        ns = 0
        esfrac = 0
        trueesfrac = 0

        for (gene,pp) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):
            n = n + 1
            ns = ns + gene_nsite[gene]
            # total post prob for positive selection thus far: estimated number of true positives in selected set
            totpp = totpp + pp
            # estimated rate of false discovery: 1 - FP / N
            efdr = 1 - totpp/n
            # estimated sensitivity (true positive rate): TP / TOT NUMBER OF POS SEL GENES
            etpr = totpp / etotp
            # estimated false negative rate
            efnr = 1 - etpr

            # estimated cumul and frac of proportion of total mass of effect size explained thus far
            cumules = cumules + gene_nsite[gene]*meanes[gene]
            esfrac = cumules / totes

            fdr = 0
            fnr = 0

            if fromsimu:
                if truees[gene] > 0:
                    tp = tp + 1
                else:
                    fp = fp + 1

                # dwestimated cumul and frac of proportion of total mass of effect size explained thus far
                truecumules = truecumules + gene_nsite[gene]*truees[gene]
                trueesfrac = truecumules / truetotes
                # false discovery rate
                fdr = fp / n
                # sensitivity or recall
                tpr = tp / ntrue
                # false negative rate
                fnr = 1 - tpr
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\t{5:7.5f}\t{6:7.5f}\t{7:7.5f}\n".format(gene,score[gene],efdr,fdr,efnr,fnr,esfrac,trueesfrac))
            else:
                outfile.write("{0:15s}\t{1:7.5f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\n".format(gene,score[gene],efdr,efnr,esfrac))

            for cutoff in cutoff_list:
                if efdr < cutoff:
                    gene_ndisc[cutoff] = n
                    gene_fdr[cutoff] = fdr
                    gene_efnr[cutoff] = efnr
                    gene_fnr[cutoff] = fnr
                    gene_cumules[cutoff] = cumules / ns * 100
                    gene_truecumules[cutoff] = truecumules / ns * 100
                    gene_esfrac[cutoff] = esfrac
                    gene_trueesfrac[cutoff] = trueesfrac

    return [gene_ndisc, gene_fdr, gene_efnr, gene_fnr, gene_cumules, gene_truecumules, gene_esfrac, gene_trueesfrac]


# collect results across all methods for a given dataset
# returns a tuple of dictionaries
# [gene_ndisc, gene_fdr, gene_efnr, gene_fnr]

def method_gene_fdr(cutoff_list, namelist, score, meanes, truees, gene_nsite, outname):

    fromsimu = len(truees)

    gene_ndisc = dict()
    gene_fdr = dict()
    gene_efnr = dict()
    gene_fnr = dict()
    gene_cumules = dict()
    gene_esfrac = dict()
    gene_truecumules = dict()
    gene_trueesfrac = dict()

    for name in namelist:
        if name == "codeml":
            for mode in ["df2", "df1", "mixdf1"]:
                [gene_ndisc[mode + "_codeml"], gene_fdr[mode + "_codeml"], gene_fnr[mode + "_codeml"], gene_truecumules[mode + "_codeml"], gene_trueesfrac[mode + "_codeml"]] = gene_codeml_fdr(cutoff_list, score["codeml"], truees, gene_nsite, outname + "_" + mode, chi_mode = mode);

        else:
            [gene_ndisc[name], gene_fdr[name], gene_efnr[name], gene_fnr[name], gene_cumules[name], gene_truecumules[name], gene_esfrac[name], gene_trueesfrac[name]] = gene_bayes_fdr(cutoff_list, score[name], meanes[name], truees, gene_nsite, outname + "_" + name);


    return {"n": gene_ndisc, 
            "fdr" : gene_fdr, 
            "efnr" : gene_efnr,
            "fnr" : gene_fnr,
            "ees" : gene_cumules,
            "es" : gene_truecumules,
            "e%es" : gene_esfrac,
            "%es" : gene_trueesfrac}


