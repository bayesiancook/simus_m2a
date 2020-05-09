import sys
import os
from numpy import mean
from numpy import product

def parse(chain_name, burnin, path = "", min_omega = 1.0) :

    current_dir = os.getcwd() + "/"
    if path != "":
        os.chdir(path)

    # open chainfile
    # get posterior mean posw and posom
    # credible interval for posom
    # keep sample of dposom
    with open(chain_name + ".chain", 'r') as chain_file:
        header = chain_file.readline()
        for i in range(burnin):
            line = chain_file.readline()

        mcmc = [line.rstrip('\n').split()[0:4] for line in chain_file]
        purom_sample = [float(sample[0]) for sample in mcmc]
        dposom_sample = [float(sample[1]) for sample in mcmc]
        purw_sample = [float(sample[2]) for sample in mcmc]
        posw_sample = [float(sample[3]) for sample in mcmc]
        sorted_dposom_sample = sorted(dposom_sample)

        # meanom = mean([ (1-pos)*(purw_sample[i]*purom_sample[i] + (1-purw_sample[i])) + pos*(1 + dposom_sample[i]) for i,pos in enumerate(posw_sample)])
        meanoma = mean([pos*dposom_sample[i] for i,pos in enumerate(posw_sample)])

        meanposw = mean(posw_sample)
        postselprob = mean([val > 0 for val in posw_sample])

        size = len(dposom_sample)
        minindex = int(0.025 * size)
        maxindex = int(0.975 * size)
        meanposom = 1 + mean(dposom_sample)
        minposom = 1 + sorted_dposom_sample[minindex]
        maxposom = 1 + sorted_dposom_sample[maxindex]


    with open(chain_name + ".param", 'r') as param_file:
        header = param_file.readline();
        [data_path, ali_file, tree_file, pi] = param_file.readline().rstrip('\n').split()

    # check for ali file and correct number of sites
    with open(data_path + ali_file, 'r') as ali_file:
        nsite = int(ali_file.readline().rstrip('\n').split()[1]) // 3

    with open(chain_name + ".sitepp", 'r') as sitepp_file:

        for i in range(burnin):
            line = sitepp_file.readline()

        mcmc = list()
        for line in sitepp_file:
            sample = line.rstrip('\n').split()
            if len(sample) != nsite :
                print("error: non matching number of sites")
                print(nsite, len(sample))
                raise
            mcmc.append(sample)
            

        if len(mcmc) != len(dposom_sample):
            print("error: mcmc samples for sitepp and gene params are not the same")
            raise

        sitepp = [mean([float(sample[i]) for sample in mcmc]) for i in range(nsite)]

        sel = dict()
        for i,pp in enumerate(sitepp):
            if pp > 0.5:
                sel[i] = pp

        postselprob2 = mean([ (dposom_sample[j] > (min_omega-1.0)) * (1 - product([1-float(sample[i]) for i in range(nsite)])) for j,sample in enumerate(mcmc)])
        postselprob3 = mean([dposom > (min_omega-1.0) for dposom in dposom_sample]) * (1 - product([1-sitepp[i] for i in range(nsite)]))

    if path != "":
        os.chdir(current_dir)

    return [postselprob, meanposw, meanposom, minposom, maxposom, sel, sitepp, postselprob2, postselprob3, meanom, meanoma]

def parse_list(basename, gene_list, burnin, path = "", min_omega = 1):

    score = dict()
    score2 = dict()
    score3 = dict()
    posw = dict()
    posom = dict()
    minposom = dict()
    maxposom = dict()
    selectedsites = dict()
    sitepp = dict()
    meanom = dict()
    meanoma = dict()

    ngene = len(gene_list)
    nsite = dict()

    for gene in gene_list:
        res = parse(basename + gene, burnin, path=path, min_omega = min_omega)
        [score[gene], posw[gene], posom[gene], minposom[gene], maxposom[gene], selectedsites[gene], sitepp[gene], score2[gene], score3[gene], meanom[gene], meanoma[gene]] = res
        nsite[gene] = len(sitepp[gene])

    totnsite = sum([ns for (gene,ns) in nsite.items()])
    grandtotom = sum([om*nsite[gene] for (gene,om) in meanom.items()])
    grandtotoma = sum([oma*nsite[gene] for (gene,oma) in meanoma.items()])
    grandmeanom = grandtotom / totnsite
    grandmeanoma = grandtotoma / totnsite
    grandmeanalpha = grandmeanoma / grandmeanom

    with open(basename + ".meanoma", 'w') as os:
        os.write("omega_tot : {0}\n".format(grandmeanom))
        os.write("omega_a   : {0}\n".format(grandmeanoma))
        os.write("alpha     : {0}\n".format(grandmeanalpha))

    with open(basename + ".sorted_oma", 'w') as outfile:

        tota = 0
        n = 0
        outfile.write("{0:18s}  {1:5s}  {2:5s}  {3:5s}  {4:5s}  {5:5s}  {6:5s}  {7:5s}  {8:5s}  {9:5s}\n".format("gene", "nsite", "om_tot", "om_a", "alpha", "%oma", "%genes", "pp", "p+", "om+"))

        for (gene,oma) in sorted(meanoma.items(), key=lambda kv: kv[1], reverse=True):
            a = oma * nsite[gene]
            tota = tota + a
            n = n + 1
            om = meanom[gene]
            alpha = oma / om
            outfile.write("{0:18s}  {1:5d}  {2:5.3f}  {3:5.3f}  {4:5.3f}  {5:5.2f}  {6:5.2f}  {7:5.2f}  {8:5.3f}  {9:5.3f}\n".format(
                gene, nsite[gene], om, oma, alpha, tota/grandtotoma, n/ngene, score[gene], posw[gene], posom[gene]))

    return [score, posw, posom, minposom, maxposom, selectedsites, sitepp, score2, score3, meanom, meanoma]


if __name__ == "__main__":

    import sys
    #chain_name = sys.argv[1]
    #burnin = int(sys.argv[2])
    #ret = parse(chain_name, burnin)
    #print(ret)

    basename = sys.argv[1]
    genefile = sys.argv[2]
    burnin = int(sys.argv[3])

    # get gene list
    with open(genefile, 'r') as listfile:
        genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

    parse_list(basename, genelist, burnin)

