import sys
import os
from numpy import mean
from numpy import product

def parse_list(chain_name, burnin, path = "", min_omega = 1.0) :

    current_dir = os.getcwd() + "/"
    if path != "":
        os.chdir(path)

    # gene general parameters
    with open(chain_name + ".param", 'r') as param_file:
        header = param_file.readline()
        [data_path, listname, tree_file] = param_file.readline().rstrip('\n').split()

    # get gene list
    # with open(listname, 'r') as listfile:
    with open(chain_name + ".genelist", 'r') as listfile:
        header = listfile.readline()
        gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for gene in listfile]

    # check for ali file and correct number of sites
    gene_nsite = dict()
    for gene in gene_list:
        with open(data_path + gene + ".ali", 'r') as ali_file:
            nsite = int(ali_file.readline().rstrip('\n').split()[1]) // 3
            gene_nsite[gene] = nsite

    ngene = len(gene_list)
    # print("number of genes : " , ngene)
    totnsite = sum([gene_nsite[gene] for gene in gene_nsite])

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

    gene_postselprob = dict()
    gene_meanposw = dict()
    gene_meanposom = dict()
    gene_minposom = dict()
    gene_maxposom = dict()
    gene_posom_sample = dict()

    for i,gene in enumerate(gene_list):
        gene_postselprob[gene] = mean([(float(sample[i]) > 0) for sample in posw_mcmc])
        gene_meanposw[gene] = mean([float(sample[i]) for sample in posw_mcmc])
        gene_meanposom[gene] = mean([float(sample[i]) for sample in posom_mcmc])

        gene_posom_sample[gene] = [float(sample[i]) for sample in posom_mcmc]
        sorted_posom_sample = sorted(gene_posom_sample[gene])
        size = len(sorted_posom_sample)
        minindex = int(0.025 * size)
        maxindex = int(0.975 * size)
        gene_minposom[gene] = sorted_posom_sample[minindex]
        gene_maxposom[gene] = sorted_posom_sample[maxindex]

    gene_sitepp = dict()
    gene_postselprob2 = dict()
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

    with open(chain_name + ".sitepp", 'r') as sitepp_file:

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

    gene_postselprob3 = dict()
    for gene in gene_list:
        gene_postselprob3[gene] = mean([posom > min_omega for posom in gene_posom_sample[gene]]) * (1 - product([1-gene_sitepp[gene][i] for i in range(gene_nsite[gene])]))
    
    if path != "":
        os.chdir(current_dir)

    return [gene_postselprob, gene_meanposw, gene_meanposom, gene_minposom, gene_maxposom, gene_selectedsites, gene_sitepp, gene_postselprob2, gene_postselprob3, hyperparams]

if __name__ == "__main__":

    import sys
    chain_name = sys.argv[1]
    burnin = int(sys.argv[2])

    ret = parse(chain_name, burnin)
    print(ret)

