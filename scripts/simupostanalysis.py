import sys
import numpy
import os
import parsecodeml
import parsem2a

current_dir = os.getcwd() + "/"
data_path = "thibotero10a_ppreduninfm2a/"
simu_folder = "thibotero10a_ppreduninfm2a/"
m2a_folder = simu_folder
codeml_folder = simu_folder

rep = 0

listname = "list10"

simu_basename = "ppreduninfm2a"
m2a_basename = "uninfm2a_"
burnin = 10

outname = simu_folder + "comp" + m2a_basename + simu_basename;

codeml_result = dict()
m2a_result = dict()
codeml_dlnl = dict()

gene_truesiteom = dict()
gene_trueposw = dict()
gene_trueposom = dict()

# get gene list
with open(listname, 'r') as listfile:
    genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

# parsing true site omega's and true param values
for gene in genelist:
        with open(current_dir + simu_folder + simu_basename + gene + "_{0}.truesiteom".format(rep)) as siteomfile:
            # header = siteomfile.readline()
            line = siteomfile.readline().rstrip('\n')
            truesiteom = [float(om) for om in line.split()]
            gene_truesiteom[gene] = truesiteom

        with open(current_dir + simu_folder + simu_basename + gene + "_{0}.trueparam".format(rep)) as trueparamfile:
            header = trueparamfile.readline()
            [true_purom, true_dposom, true_purw, true_posw] = [float(item) for item in trueparamfile.readline().rstrip('\n').split()[0:4]]
            gene_trueposw[gene] = true_posw
            gene_trueposom[gene] = 1.0 + true_dposom

# parsing codeml results
for gene in genelist:
    codeml_res = parsecodeml.parse(current_dir + codeml_folder + simu_basename + gene + "_{0}.codeml".format(rep))
    codeml_result[gene] = codeml_res
    codeml_dlnl[gene] = codeml_res[1] - codeml_res[0]

# parsing m2a results
os.chdir(current_dir + m2a_folder)
for gene in genelist:
    m2a_res = parsem2a.parse(m2a_basename + simu_basename + gene + "_{0}".format(rep), burnin)
    m2a_result[gene] = m2a_res
os.chdir("..")

cutoff_list = [0.5, 0.6, 0.7, 0.8, 0.9]

totbayes_ndisc = dict()
totbayes_fdr = dict()
totbayes_efdr = dict()
totml_ndisc = dict()
totml_fdr = dict()
totml_efdr = dict()

for cutoff in cutoff_list:
    totbayes_ndisc[cutoff] = 0
    totbayes_fdr[cutoff] = 0
    totbayes_efdr[cutoff] = 0
    totml_ndisc[cutoff] = 0
    totml_fdr[cutoff] = 0
    totml_efdr[cutoff] = 0

with open(outname + ".sortedgenes", 'w') as geneoutfile:
    with open(outname + ".sitepp", 'w') as siteoutfile:

        for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):

            bayes_posw = m2a_result[gene][1]
            bayes_posom = m2a_result[gene][2]
            (bayes_min, bayes_max) = m2a_result[gene][3]

            ml_posw = codeml_result[gene][6]
            ml_posom = codeml_result[gene][9]

            bayes_sitepp = m2a_result[gene][4]
            bayes_siteom = m2a_result[gene][5]
            codeml_sites = codeml_result[gene][11]

            nsite = len(bayes_sitepp)
            ml_sitepp = [0 for i in range(nsite)]
            ml_siteom = [0 for i in range(nsite)]
            # comp_sites = dict()
            for i in codeml_sites:
                (ml_pp, ml_om) = codeml_sites[i]
                # bayes_pp = bayes_sitepp[i-1]
                # bayes_om = bayes_siteom[i-1]
                ml_sitepp[i-1] = float(ml_pp)
                ml_siteom[i-1] = float(ml_pp)
                # comp_sites[i] = (int(100*ml_pp), int(100*bayes_pp))
                # siteoutfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(gene,ml_pp,bayes_pp,ml_om,bayes_om))

            geneoutfile.write("{0}\t{1:6.2f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.2f}\t{6:6.2f}\t{7:6.2f} ({8:6.2f},{9:6.2f})\t".format(gene, lnl, gene_trueposw[gene], ml_posw, bayes_posw, gene_trueposom[gene], ml_posom, bayes_posom, bayes_min, bayes_max))

            for cutoff in cutoff_list:
                bayes_ndisc = sum([(pp > cutoff) for pp in bayes_sitepp])
                bayes_efdr = sum([pp for pp in bayes_sitepp if pp > cutoff])
                bayes_fdr = sum([ ((pp > cutoff) and (gene_truesiteom[gene][i] > 1.0)) for i,pp in enumerate(bayes_sitepp) ])

                totbayes_ndisc[cutoff] = totbayes_ndisc[cutoff] + bayes_ndisc
                totbayes_efdr[cutoff] = totbayes_efdr[cutoff] + bayes_efdr
                totbayes_fdr[cutoff] = totbayes_fdr[cutoff] + bayes_fdr

                ml_ndisc = sum([(pp > cutoff) for pp in ml_sitepp])
                ml_efdr = sum([pp for pp in ml_sitepp if pp > cutoff])
                ml_fdr = sum([ ((pp > cutoff) and (gene_truesiteom[gene][i] > 1.0)) for i,pp in enumerate(ml_sitepp) ])

                totml_ndisc[cutoff] = totml_ndisc[cutoff] + ml_ndisc
                totml_efdr[cutoff] = totml_efdr[cutoff] + ml_efdr
                totml_fdr[cutoff] = totml_fdr[cutoff] + ml_fdr

                geneoutfile.write("{0:5d} {1:5d} {2:6.4f}".format(bayes_ndisc, bayes_fdr, bayes_efdr))
                geneoutfile.write("{0:5d} {1:5d} {2:6.4f}".format(ml_ndisc, ml_fdr, ml_efdr))

            geneoutfile.write("\n")

print()
print("cutoff\tbayes\tefdr\tfdr\tml\tefdr\tfdr")
for cutoff in cutoff_list:

    bayes_fdr = 0
    bayes_efdr = 0
    if totbayes_ndisc[cutoff]:
        bayes_fdr = totbayes_fdr[cutoff] / totbayes_ndisc[cutoff]
        bayes_efdr = totbayes_efdr[cutoff] / totbayes_ndisc[cutoff]

    ml_fdr = 0
    ml_efdr = 0
    if totml_ndisc[cutoff]:
        ml_fdr = totml_fdr[cutoff] / totml_ndisc[cutoff]
        ml_efdr = totml_efdr[cutoff] / totml_ndisc[cutoff]

    print("{0:5.2f}\t{1:5d}\t{2:5.2f}\t{3:5.2f}\t{4:5d}\t{5:5.2f}\t{6:5.2f}".format(cutoff, totbayes_ndisc[cutoff], bayes_efdr, bayes_fdr, totml_ndisc[cutoff], ml_efdr, ml_fdr))

print()



                
