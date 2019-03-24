#! /usr/bin/python3.5

import sys
import numpy
import os
import datalib

def compare(run_names, gene_list_name, taxlist1, taxlist2, outname, data_path = "", target_path = ""):

    data_dir = data_path + "/"
    target_dir = target_path + "/"

    # get gene list
    with open(gene_list_name, 'r') as list_file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in list_file]
    ngene = len(gene_list)

    # check for ali file and correct number of sites
    gene_nsite = dict()
    for gene in gene_list:
        with open(data_dir + gene + ".ali", 'r') as ali_file:
            nsite = int(ali_file.readline().rstrip('\n').split()[1]) // 3
            gene_nsite[gene] = nsite

    # get the two taxon lists
    with open(taxlist1, 'r') as taxfile:
        taxon_list1 = [line.rstrip('\n') for line in taxfile]
    with open(taxlist2, 'r') as taxfile:
        taxon_list2 = [line.rstrip('\n') for line in taxfile]

    # initialize dictionaries
    all_pp = dict()
    all_posw = dict()
    all_posom = dict()
    all_sitepp = dict()

    # read all gene and site post probs
    for name in run_names:
        print(name)

        gene_pp = dict()
        gene_posw = dict()
        gene_posom = dict()

        with open(name + ".genepost", 'r') as infile:
            header = infile.readline()
            for line in infile:
                [gene, pp, posw, posom] = line.rstrip('\n').split()
                gene_pp[gene] = float(pp)
                gene_posw[gene] = float(posw)
                gene_posom[gene] = float(posom)

        all_pp[name] = gene_pp
        all_posw[name] = gene_posw
        all_posom[name] = gene_posom

        gene_sitepp = dict()
        for gene in gene_list:
            gene_sitepp[gene] = [0 for site in range(gene_nsite[gene])]

        with open(name + ".sitepost", 'r') as infile:
            header = infile.readline()
            for line in infile:
                [gene, site, pp] = line.rstrip('\n').split()
                gene_sitepp[gene][int(site)] = float(pp)

        all_sitepp[name] = gene_sitepp
            
    with open(target_dir + outname + ".comp", 'w') as outfile:

        c = 0.95

        # compare gene post probs
        for i,a in enumerate(run_names):
            for j,b in enumerate(run_names):

                if i != j :
                    comp = 0
                    for gene in gene_list:
                        if (all_pp[a][gene] > c) and (all_pp[b][gene] > c):
                            comp += 1
                    outfile.write("{0:6d}\t".format(comp))

                else :
                    comp = 0
                    for gene in gene_list:
                        if all_pp[a][gene] > c:
                            comp += 1
                    outfile.write("{0:6d}\t".format(comp))
            outfile.write("\n")

        outfile.write("\n")
        # compare site post probs
        # marginally
        for i,a in enumerate(run_names):
            for j,b in enumerate(run_names):

                if i != j :
                    comp = 0
                    for gene in gene_list:
                        for site in range(gene_nsite[gene]):
                            if (all_sitepp[a][gene][site] > c) and (all_sitepp[b][gene][site] > c):
                                comp += 1
                    outfile.write("{0:6d}\t".format(comp))

                else :
                    comp = 0
                    for gene in gene_list:
                        for site in range(gene_nsite[gene]):
                            if all_sitepp[a][gene][site] > c:
                                comp += 1
                    outfile.write("{0:6d}\t".format(comp))
            outfile.write("\n")

        alphabet = ['AAC', 'TGG', 'CGA', '---', '---', '---', '---']
        # output alignments for genes that are congruent
        gene_comp = dict()
        for gene in gene_list:
            congruent = numpy.prod([all_pp[name][gene] > c for name in run_names])
            if congruent:
                sitepp = [numpy.sum([all_sitepp[name][gene][site] > c for name in run_names]) for site in range(gene_nsite[gene])]
                n = sum([pp>1 for pp in sitepp])
                if n > 3:
                    ppseq = ''.join([alphabet[x] for x in sitepp])
                    (ntaxa, nsite, ali) = datalib.read_phylip(data_dir + gene + ".ali")
                    with open(target_dir + "pos_" + gene + ".fas", 'w') as outfile:
                        outfile.write(">PP_CONS\n{0}\n".format(ppseq))
                        for tax in taxon_list1:
                            outfile.write(">{0}\n{1}\n".format(tax,ali[tax]))
                        for tax in taxon_list2:
                            outfile.write(">{0}\n{1}\n".format(tax,ali[tax]))
                        # for (tax,seq) in ali.items():
                        #    outfile.write(">{0}\n{1}\n".format(tax,seq))

# main
dir_list = ["multiommv10_euarch", "multiommv10_laura"]
# dir_list = ["multiommv10_placnr", "multiommv10_euarch", "multiommv10_laura"]
basename_list = ["indmm2a"]
run_names = [ dir + "/" + name for dir in dir_list for name in basename_list]
compare(run_names, "all.list", "euarch.taxlist", "laura.taxlist", "compruns", data_path = "multiommv10_placnr", target_path = "results")

