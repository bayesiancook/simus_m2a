#! /usr/bin/python3.5

import sys
import os
import numpy

def get_true_params(exp_folder, min_omega = 1):

    exp_dir = exp_folder + "/"
    listname = "all.list"

    # get gene list
    with open(exp_dir + listname, 'r') as listfile:
        genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

    ngene = len(genelist)

    # get number of taxa and sites per gene
    gene_ntaxa = dict()
    gene_nsite = dict()
    for gene in genelist:
        with open(exp_dir + "data/" + gene + ".ali") as genefile:
            [ntaxa, nsite] = genefile.readline().rstrip('\n').split()
            gene_ntaxa[gene] = int(ntaxa)
            gene_nsite[gene] = int(nsite) // 3

    truepurw = dict()
    trueposw = dict()
    truepos2 = dict()
    truepurom = dict()
    truedposom = dict()
    truesiteom = dict()
    for gene in genelist:
        with open(exp_dir + "data/" + gene + ".trueparam") as genefile:
            header = genefile.readline()
            [purom, dposom, purw, posw] = [float(f) for f in genefile.readline().rstrip('\n').split()[0:4]]
            truepurom[gene] = purom
            truedposom[gene] = dposom
            truepurw[gene] = purw
            trueposw[gene] = posw
        with open(exp_dir + "data/" + gene + ".truesiteom") as genefile:
            truesiteom[gene] = [float(f) for f in genefile.readline().rstrip('\n').split()]
            if len(truesiteom[gene]) != gene_nsite[gene]:
                print("error: non matching number of sites :", gene_nsite[gene], len(truesiteom[gene]))
                sys.exit()
        truepos2[gene] = sum(om > min_omega for om in truesiteom[gene])

    return [truepurw, trueposw, truepurom, truedposom, truesiteom, truepos2]

def get_empirical_hyperparams(truepurw, trueposw, truepurom, truedposom):

    purw_mean = numpy.mean([x for (gene,x) in truepurw.items()])
    purw_var = numpy.var([x for (gene,x) in truepurw.items()])
    purw_invconc = 1.0 / (purw_mean * (1-purw_mean) / purw_var - 1)

    posw_mean = numpy.mean([x for (gene,x) in trueposw.items() if x > 0])
    posw_var = numpy.var([x for (gene,x) in trueposw.items() if x > 0])
    posw_invconc = 1.0 / (posw_mean * (1-posw_mean) / posw_var - 1)

    purom_mean = numpy.mean([x for (gene,x) in truepurom.items()])
    purom_var = numpy.var([x for (gene,x) in truepurom.items()])
    purom_invconc = 1.0 / (purom_mean * (1-purom_mean) / purom_var - 1)

    dposom_mean = numpy.mean([x for (gene,x) in truedposom.items() if trueposw[gene] > 0])
    dposom_var = numpy.var([x for (gene,x) in truedposom.items() if trueposw[gene] > 0])
    dposom_invshape = dposom_var / dposom_mean / dposom_mean

    pi = numpy.mean([x>0 for (gene,x) in trueposw.items()])

    return (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape)

def get_empirical_moments(exp_folder):

    [truepurw, trueposw, truepurom, truedposom, truesiteom, truepos2] = get_true_params(exp_folder)
    return get_empirical_hyperparams(truepurw, trueposw, truepurom, truedposom)

if __name__ == "__main__":

    if len(sys.argv) == 1:
        print("simuparams.py experiment")
        sys.exit()

    exp_folder = sys.argv[1]

    [truepurw, trueposw, truepurom, truedposom, truesiteom, truepos2] = get_true_params(exp_folder)

    (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = get_empirical_hyperparams(truepurw, trueposw, truepurom, truedposom)

    print("{0:15s}\t{1:>5s}\t{2:>5s}\t{3:>5s}\t{4:>5s}\t{5:>5s}".format("gene","purw","posw","purom","dpsom", "pos2"))

    for (gene,posw) in sorted(trueposw.items(), key=lambda kv: kv[1], reverse=True):
        print("{0}\t{1:5.3f}\t{2:5.3f}\t{3:5.2f}\t{4:5.2f}\t{5:5d}".format(gene, truepurw[gene], trueposw[gene], truepurom[gene], 1.0 + truedposom[gene], truepos2[gene]))

    print()
    print("summary across genes:")
    print()
    print("pi     : {0:5.3f}".format(pi))
    print("purw   : {0:5.3f}\t{1:5.3f}\t{2:5.3f}".format(purw_mean, numpy.sqrt(purw_var), purw_invconc))
    print("posw   : {0:5.3f}\t{1:5.3f}\t{2:5.3f}".format(posw_mean, numpy.sqrt(posw_var), posw_invconc))
    print("purom  : {0:5.3f}\t{1:5.3f}\t{2:5.3f}".format(purom_mean, numpy.sqrt(purom_var), purom_invconc))
    print("dposom : {0:5.3f}\t{1:5.3f}\t{2:5.3f}".format(dposom_mean, numpy.sqrt(dposom_var), dposom_invshape))
    print("pos2   : {0:5.3f}".format(sum([truepos2[gene]>0 for gene in truepos2]) / len(truepos2)))
    print()
 

