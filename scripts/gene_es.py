#! /usr/bin/python3.5

import sys
import numpy
import os

default_cutoff_list = [0.1, 0.2, 0.3]

def gene_es(score, meanes, truees, nsite, outname, cutoff_list = default_cutoff_list):

    gene_ndisc = dict()
    gene_cumules = dict()
    gene_truecumules = dict()
    gene_trueesfrac = dict()
    gene_efdr = dict()
    gene_fdr = dict()
    for cutoff in cutoff_list:
        gene_ndisc[cutoff] = 0
        gene_cumules[cutoff] = 0
        gene_truecumules[cutoff] = 0
        gene_trueesfrac[cutoff] = 0
        gene_efdr[cutoff] = 0
        gene_fdr[cutoff] = 0

    ngene = len(meanes)
    fromsimu = len(truees)

    with open(outname + ".sortedes", 'w') as outfile:

        if fromsimu:

            outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}  {7:5s}  {8:5s}  {9:5s}\n".format("#gene", "nsite", "ees", "es", "e%es", "%es", "%genes", "pp", "efdr", "fdr"))

            totnsite = sum([ns for (gene,ns) in nsite.items()])
            totes = sum([es*nsite[gene] for (gene,es) in meanes.items()])
            tottruees = sum([es*nsite[gene] for (gene,es) in truees.items()])
            truecumules = 0
            cumules = 0
            n = 0
            ns = 0
            fp = 0
            tp = 0
            totpp = 0

            for (gene,es) in sorted(meanes.items(), key=lambda kv: kv[1], reverse=True):

                n = n + 1
                ns = ns + nsite[gene]

                totpp = totpp + score[gene]
                # estimated rate of false discovery: 1 - FP / N
                efdr = 1 - totpp/n

                if truees[gene] > 0:
                    tp = tp + 1
                else:
                    fp = fp + 1
                fdr = fp / n

                cumules = cumules + es*nsite[gene]
                esfrac = cumules / totes
                truecumules = truecumules + truees[gene]*nsite[gene]
                trueesfrac = truecumules / tottruees

                outfile.write("{0:18s}  {1:5d}  {2:5.3f}  {3:5.3f}  {4:5.2f}  {5:5.2f}  {6:5.2f}  {7:5.3f}  {8:5.3f}  {9:5.3f}\n".format(
                    gene, nsite[gene], es, truees[gene], esfrac, trueesfrac, n/ngene, score[gene], efdr, fdr))

                for cutoff in cutoff_list:
                    if esfrac < 1-cutoff:
                        gene_ndisc[cutoff] = n
                        gene_cumules[cutoff] = cumules / ns * 100
                        gene_truecumules[cutoff] = truecumules / ns * 100
                        gene_trueesfrac[cutoff] = trueesfrac
                        gene_efdr[cutoff] = efdr
                        gene_fdr[cutoff] = fdr

        else:

            outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}\n".format("#gene", "nsite", "ees", "e%es", "%genes", "pp", "efdr"))

            totnsite = sum([ns for (gene,ns) in nsite.items()])
            totes = sum([es*nsite[gene] for (gene,es) in meanes.items()])
            cumules = 0
            n = 0
            ns = 0
            totpp = 0

            for (gene,es) in sorted(meanes.items(), key=lambda kv: kv[1], reverse=True):

                n = n + 1
                ns = ns + nsite[gene]

                totpp = totpp + score[gene]
                # estimated rate of false discovery: 1 - FP / N
                efdr = 1 - totpp/n

                cumules = cumules + es*nsite[gene]
                esfrac = cumules / totes

                outfile.write("{0:18s}  {1:5d}  {2:5.3f}  {3:5.3f}  {4:5.2f}  {5:5.2f}  {6:5.2f}\n".format(
                    gene, nsite[gene], es, esfrac, n/ngene, score[gene], efdr))

                for cutoff in cutoff_list:
                    if esfrac < 1-cutoff:
                        gene_ndisc[cutoff] = n
                        gene_cumules[cutoff] = cumules / ns * 100
                        gene_efdr[cutoff] = efdr

    return [gene_ndisc, gene_cumules, gene_truecumules, gene_trueesfrac, gene_efdr, gene_fdr]



