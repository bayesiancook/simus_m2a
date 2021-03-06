#! /beegfs/data/soft/anaconda3/bin/python3.7

import sys
import os
from numpy import mean
from numpy import product
from fdr import gene_bayes_fdr
from gene_es import gene_es

hypernamelist = ['purom_mean', 'purom_invconc', 'dposom_mean', 'dposom_invshape', 'purw_mean', 'purw_invconc', 'posw_mean', 'posw_invconc', 'pi']

def parse_list(chain_name, burnin, with_sites = True, write_output = False, path = "", min_omega = 1.0) :

    current_dir = os.getcwd() + "/"
    if path != "":
        os.chdir(path)

    print("get params and gene lists")
    # gene general parameters
    with open(chain_name + ".param", 'r') as param_file:
        header = param_file.readline()
        [data_path, listname, tree_file] = param_file.readline().rstrip('\n').split()

    # get gene list
    with open(chain_name + ".genelist", 'r') as listfile:
        header = listfile.readline()
        ngene = int(header.rstrip('\n').split()[0])
        gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for i,gene in enumerate(listfile) if i < ngene]

    gene_nsite = dict()

    # get original gene list
    with open(listname, 'r') as listfile:
        header = listfile.readline().rstrip('\n')
        if header == "ALI":
            header = listfile.readline().rstrip('\n')
            ngene = int(header.split()[0])
            original_gene_list = []
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
            for gene in gene_list:
                with open(data_path + gene + ".ali", 'r') as ali_file:
                    nsite = int(ali_file.readline().rstrip('\n').split()[1]) // 3
                    gene_nsite[gene] = nsite

    ngene = len(gene_list)
    # print("number of genes : " , ngene)
    totnsite = sum([gene_nsite[gene] for gene in gene_nsite])

    print("processing hyperparams")
    # open chain file and get hyperparams
    with open(chain_name + ".chain", 'r') as chain_file:
        # header = chain_file.readline()
        for i in range(burnin):
            line = chain_file.readline()

        mcmc = [line.rstrip('\n').split()[0:9] for line in chain_file]

        hyperparams = list()
        for i in range(9):
            sorted_sample = sorted([float(s[i]) for s in mcmc])
            post_mean = mean(sorted_sample)
            size = len(sorted_sample)
            minindex = int(0.025 * size)
            maxindex = int(0.975 * size)
            post_min = sorted_sample[minindex]
            post_max = sorted_sample[maxindex]
            hyperparams.append((post_mean, post_min, post_max))

    print("processing gene pps")
    # open posom and posw files 
    with open(chain_name + ".posw", 'r') as posw_file:
        # header = posw_file.readline()
        for i in range(burnin):
            line = posw_file.readline()

        posw_mcmc = [line.rstrip('\n').split()[0:ngene] for line in posw_file]

    with open(chain_name + ".posom", 'r') as posom_file:
        # header = posom_file.readline()

        for i in range(burnin):
            line = posom_file.readline()

        posom_mcmc = [line.rstrip('\n').split()[0:ngene] for line in posom_file]

    print("post processing gene pps")
    gene_postselprob = dict()
    gene_meanposw = dict()
    gene_meanposom = dict()
    gene_meanoma = dict()
    gene_minposom = dict()
    gene_maxposom = dict()
    gene_posom_sample = dict()

    for i,gene in enumerate(gene_list):
        gene_postselprob[gene] = mean([(float(sample[i]) > 0) for sample in posw_mcmc])
        gene_meanposw[gene] = mean([float(sample[i]) for sample in posw_mcmc])
        gene_meanposom[gene] = mean([float(sample[i]) for sample in posom_mcmc])
        gene_meanoma[gene] = mean([float(posw[i]) * (float(posom_mcmc[j][i]) - 1) for j,posw in enumerate(posw_mcmc)])

        gene_posom_sample[gene] = [float(sample[i]) for sample in posom_mcmc]
        sorted_posom_sample = sorted(gene_posom_sample[gene])
        size = len(sorted_posom_sample)
        minindex = int(0.025 * size)
        maxindex = int(0.975 * size)
        gene_minposom[gene] = sorted_posom_sample[minindex]
        gene_maxposom[gene] = sorted_posom_sample[maxindex]

    gene_sitepp = dict()
    gene_postselprob2 = dict()
    gene_postselprob3 = dict()

    for gene in gene_list:
        gene_sitepp[gene] = [0 for site in range(gene_nsite[gene])]
        gene_postselprob2[gene] = 0

    gene_selectedsites = dict()

    nsample = 0
    for gene in gene_list:
        if nsample == 0:
            nsample = len(gene_posom_sample[gene])
        else:
            if nsample != len(gene_posom_sample[gene]):
                print("error: non matching sample size")
                raise

    if with_sites:
        with open(chain_name + ".sitepp", 'r') as sitepp_file:

            print("processing sitepp file")
            for i in range(burnin):
                line = sitepp_file.readline()

            sample_size = 0
            for line in sitepp_file:
                sample = line.rstrip('\n').split()
                if len(sample) != totnsite + ngene:
                    print("error: non matching number of sites")
                    print(totnsite, ngene, len(sample))
                    raise
                
                index = 0
                # for j in range(ngene):
                for j,gene in enumerate(gene_list):
                    name = sample[index].replace(".ali","")
                    if name != gene:
                        print("error: non matching gene name")
                        sys.exit()

                    index = index + 1

                    sitepp = gene_sitepp[gene]

                    pp_prod = 1.0
                    for k in range(gene_nsite[gene]):
                        pp = float(sample[index+k])
                        sitepp[k] = sitepp[k] + pp
                        pp_prod = pp_prod * (1-pp)
                    index = index + gene_nsite[gene]

                    gene_postselprob2[gene] = gene_postselprob2[gene] + (1 - pp_prod) * (gene_posom_sample[gene][sample_size] > min_omega)

                sample_size = sample_size + 1

            print("post processing site pps")
            for gene in gene_list:
                gene_postselprob2[gene] = gene_postselprob2[gene] / sample_size
                sitepp = gene_sitepp[gene]
                for k in range(gene_nsite[gene]):
                    sitepp[k] = sitepp[k] / sample_size

                sel = dict()
                for i,pp in enumerate(sitepp):
                    if pp > 0.5:
                        sel[i] = pp

                gene_selectedsites[gene] = sel

        for gene in gene_list:
            gene_postselprob3[gene] = mean([posom > min_omega for posom in gene_posom_sample[gene]]) * (1 - product([1-gene_sitepp[gene][i] for i in range(gene_nsite[gene])]))
    
    if write_output:
        with open(chain_name + ".postanalysis", 'w') as outfile:
            outfile.write("gene\tpp\tposw\tom")
            if with_sites:
                outfile.write("\tsites")
            outfile.write("\n")
            for gene in original_gene_list:
                outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(gene, gene_postselprob[gene], gene_meanposw[gene], gene_meanposom[gene], gene_minposom[gene], gene_maxposom[gene]))
                if with_sites:
                    for k in range(gene_nsite[gene]):
                        outfile.write("\t{0}".format(int(100 * gene_sitepp[gene][k])))
                outfile.write("\n")


        totnsite = sum([ns for (gene,ns) in gene_nsite.items()])
        grandtotoma = sum([oma*gene_nsite[gene] for (gene,oma) in gene_meanoma.items()])
        grandmeanoma = grandtotoma / totnsite

        with open(chain_name + ".sortedoma", 'w') as outfile:

            tota = 0
            n = 0
            outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}  {7:5s}\n".format("gene", "nsite", "om_a", "%oma", "%genes", "pp", "p+", "om+"))

            for (gene,oma) in sorted(gene_meanoma.items(), key=lambda kv: kv[1], reverse=True):
                a = oma * gene_nsite[gene]
                tota = tota + a
                n = n + 1
                outfile.write("{0:18s}  {1:5d}  {2:5.3f}  {3:5.2f}  {4:5.2f}  {5:5.3f}  {6:5.3f}  {7:5.3f}\n".format(
                    gene, gene_nsite[gene], oma, tota/grandtotoma, n/ngene, gene_postselprob[gene], gene_meanposw[gene], gene_meanposom[gene]))

    if path != "":
        os.chdir(current_dir)

    return [gene_postselprob, gene_meanposw, gene_meanposom, gene_minposom, gene_maxposom, gene_selectedsites, gene_sitepp, gene_postselprob2, gene_postselprob3, hyperparams, gene_meanoma, gene_nsite]

if __name__ == "__main__":

    import sys
    if len(sys.argv) == 1:
        print("parsemm2a chain_name burnin [-s]")
        sys.exit()

    chain_name = sys.argv[1]
    burnin = int(sys.argv[2])
    with_sites = True
    if sys.argv[3] == "-s":
        with_sites = False

    res = parse_list(chain_name, burnin, with_sites=with_sites, write_output = True)
    [score, posw, posom, minposom, maxposom, selectedsites, sitepp, score2, score3, hyperparams, meanoma, gene_nsite] = res[0:12]

    gene_bayes_fdr([0.05, 0.1, 0.3], score, posw, dict(), posw, dict(), gene_nsite, chain_name)
    print("fdr series in " + chain_name + ".genefdr")


    print("oma series in " + chain_name + ".geneoma")
    with open(chain_name+ ".posthyper", 'w') as outfile:
        for i,name in enumerate(hypernamelist):
            outfile.write("{0:20s}".format(name))
            (mean,min,max) = hyperparams[i]
            outfile.write("\t{0:6.4f} ({1:6.4f} , {2:6.4f})\n".format(mean,min,max))

    print("estimated hyperparameters in " + chain_name + ".posthyper")

    gene_es(score, posw, dict(), gene_nsite, chain_name, cutoff_list = [0.1, 0.2, 0.3])

