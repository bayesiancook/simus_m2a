#! /usr/bin/python3.5

import sys
import os
from numpy import mean
from numpy import product
from fdr import gene_fdr

def parse_list(chain_name, burnin, write_output = False, path = "", min_omega = 1.0) :

    current_dir = os.getcwd() + "/"
    if path != "":
        os.chdir(path)

    print("get params and gene lists")
    # gene general parameters
    with open(chain_name + ".param", 'r') as param_file:
        header = param_file.readline()
        [listname, tree_file] = param_file.readline().rstrip('\n').split()

    # get gene list
    with open(chain_name + ".genelist", 'r') as listfile:
        header = listfile.readline()
        ngene = int(header.rstrip('\n').split()[0])
        gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for i,gene in enumerate(listfile) if i < ngene]

    # get original gene list
    with open(listname, 'r') as listfile:
        header = listfile.readline().rstrip('\n')
        if header == "ALI":
            header = listfile.readline().rstrip('\n')
            ngene = int(header.split()[0])
            original_gene_list = []
            gene_nsite = dict()
            for i in range(ngene):
                gene = listfile.readline().rstrip('\n').replace(".ali","")
                original_gene_list.append(gene)
                (ntax,npos) = listfile.readline().rstrip('\n').split()
                ntaxa = int(ntax)
                for j in range(ntaxa):
                    line = listfile.readline()
                    if not j:
                        (tax,seq) = line.rstrip('\n').split()
                        gene_nsite[gene] = len(seq) // 3
        else:
            ngene = int(header.split()[0])
            original_gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for i,gene in enumerate(listfile) if i < ngene]
            # original_gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for gene in listfile]

            # check for ali file and correct number of sites
            gene_nsite = dict()
            for gene in gene_list:
                with open(data_path + gene + ".ali", 'r') as ali_file:
                    nsite = int(ali_file.readline().rstrip('\n').split()[1]) // 3
                    gene_nsite[gene] = nsite

    ngene = len(gene_list)
    # print("number of genes : " , ngene)
    totnsite = sum([gene_nsite[gene] for gene in gene_nsite])

    print("processing gene oms")
    # open posom files 
    with open(chain_name + ".geneom", 'r') as posom_file:
        # header = posom_file.readline()
        for i in range(burnin):
            line = posom_file.readline()

        posom_mcmc = [line.rstrip('\n').split()[0:ngene] for line in posom_file]

    print("post processing gene oms")
    gene_postselprob = dict()
    gene_meanposom = dict()
    gene_minposom = dict()
    gene_maxposom = dict()
    gene_posom_sample = dict()

    alpha = 0.05
    for i,gene in enumerate(gene_list):
        gene_postselprob[gene] = mean([(float(sample[i]) > min_omega) for sample in posom_mcmc])
        gene_meanposom[gene] = mean([float(sample[i]) for sample in posom_mcmc])
        gene_posom_sample[gene] = [float(sample[i]) for sample in posom_mcmc]
        sorted_posom_sample = sorted(gene_posom_sample[gene])
        size = len(sorted_posom_sample)
        minindex = int(alpha / 2 * size)
        maxindex = int((1 - alpha) / 2 * size)
        gene_minposom[gene] = sorted_posom_sample[minindex]
        gene_maxposom[gene] = sorted_posom_sample[maxindex]

    if write_output:
        with open(chain_name + ".postanalysis", 'w') as outfile:
            outfile.write("gene\tpp\tom\tmin\tmax\n")
            for gene in original_gene_list:
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene, gene_postselprob[gene], gene_meanposom[gene], gene_minposom[gene], gene_maxposom[gene]))

    if path != "":
        os.chdir(current_dir)

    return [gene_postselprob, gene_meanposom, gene_minposom, gene_maxposom]

if __name__ == "__main__":

    import sys
    if len(sys.argv) == 1:
        print("parsemmutsel chain_name burnin [-s]")
        sys.exit()

    chain_name = sys.argv[1]
    burnin = int(sys.argv[2])
    min_omega = float(sys.argv[3])
    res = parse_list(chain_name, burnin, write_output = True, min_omega = min_omega)

    [score, posom, minposom, maxposom] = res[0:4]
    truepos = dict()
    cutoff_list = [0.5, 0.7, 0.9]
    [gene_ndisc, gene_fp, gene_efdr, gene_etpr] = gene_fdr(cutoff_list, score, truepos, chain_name)

    with open(chain_name + ".genefdr", 'w') as outfile:
        outfile.write("{0:5s} {1:5s} {2:5s} {3:5s}\n".format("c", "ndisc", "efdr", "etpr"))
        for cutoff in cutoff_list:
            ndisc = gene_ndisc[cutoff]
            if ndisc:
                efdr = gene_efdr[cutoff]
                etpr = gene_etpr[cutoff]
                outfile.write("{0:5.1f} {1:5d} {2:5.2f} {3:5.2f}\n".format(cutoff, ndisc, efdr, etpr))


    print("stats in " + chain_name + ".genefdr")


