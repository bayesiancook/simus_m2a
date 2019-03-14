import sys
import os
import re
import runmm2a
import datalib
from numpy import random

def single_gene_postpred(exp_folder, listname, basename, target_folder, prop_pos = 0.20, shrink_posw = 1.0, shrink_dposom = 1.0, prob_posw = 0, prob_dposom = 0, burnin = 10, every = 10):

    nrep = 1
    until = burnin + nrep * every + 1
    rep = 0

    exp_dir = exp_folder + "/"
    single_dir = exp_dir + "singlegene/"

    # this should normally be replicated
    target_dir = target_folder + "/"
    # simus in data sub-folder of target folder
    simu_folder = target_dir + "data"
    simu_dir = simu_folder + "/"

    if os.path.exists(target_folder):
        print("target folder ({0}) already exists".format(target_folder))
        raise

    os.system("mkdir " + target_folder)
    os.system("mkdir " + simu_folder)
    os.system("cp " + exp_dir + "all.tree " + target_folder)
    os.system("cp " + exp_dir + "all.list " + target_folder)

    match_gene = r"^gene\s+"

    genelist = dict()
    with open(listname) as listfile:

        for line in listfile:
            m = re.match(match_gene,line.rstrip('\n'))
            if m:
                break

        for line in listfile:
            gene2lnl = line.rstrip('\n').split()[0:2]
            genelist[gene2lnl[0]] = float(gene2lnl[1])

    ngene = len(genelist)
    ngene_pos = int(ngene * prop_pos)

    # make ppred simulations for all genes
    gene_index = 0
    for (gene,lnl) in sorted(genelist.items(), key=lambda kv: kv[1], reverse=True):

        # check trees are there
        if not os.path.exists(single_dir + gene + ".tree"):
            print("gene tree " + single_dir + gene + ".tree does not exist")
            raise
        os.system("cp " + single_dir + gene + ".tree " + simu_dir)

        # readcodonm2a with ppred option
        gene_index = gene_index + 1
        ppred_command = "readcodonm2a -x {0} {1} {2} -ppred -null {3}{4}".format(burnin,every,until,basename,gene)
        if gene_index <= ngene_pos:
            sposw = 1.0
            if random.random() < prob_posw:
                sposw = shrink_posw
            sdposom = 1.0
            if random.random() < prob_dposom:
                sdposom = shrink_dposom

            ppred_command = "readcodonm2a -x {0} {1} {2} -ppred -shrinkposw {3} -shrinkdposom {4} {5}{6}".format(burnin, every, until, sposw, sdposom, basename, gene)

        # execute ppred and move resulting files in simu folder
        cwd = os.getcwd()
        os.chdir(single_dir)
        os.system(ppred_command)
        os.chdir(cwd)
        os.system("mv {0}ppred{1}{2}_{3}.ali {4}{2}.ali".format(single_dir, basename, gene, rep, simu_dir))
        os.system("mv {0}ppred{1}{2}_{3}.trueparam {4}{2}.trueparam".format(single_dir, basename, gene, rep, simu_dir))
        os.system("mv {0}ppred{1}{2}_{3}.truesiteom {4}{2}.truesiteom".format(single_dir, basename, gene, rep, simu_dir))


# assumes post pred has already been run
def multi_gene_postpred(exp_folder, basename, target_folder):

    exp_dir = exp_folder + "/"
    multi_dir = exp_dir + "multigene/"

    # this should normally be replicated
    target_dir = target_folder + "/"
    # simus in data sub-folder of target folder
    simu_folder = target_dir + "data"
    simu_dir = simu_folder + "/"

    if os.path.exists(target_folder):
        print("target folder ({0}) already exists".format(target_folder))
        raise

    os.system("mkdir " + target_folder)
    os.system("mkdir " + simu_folder)
    os.system("cp " + exp_dir + "all.tree " + target_folder)
    os.system("cp " + exp_dir + "all.list " + target_folder)

    with open(exp_dir + "all.list", 'r') as listfile:
        gene_list = [gene.rstrip('\n').split()[0].replace(".ali","") for gene in listfile]

    for gene in gene_list:
        # open single gene alignment from original data
        (ntaxa, nsite, original_ali) = datalib.read_phylip(exp_dir + "data/" + gene + ".ali")
        #taxlist = [tax for (tax,seq) in original_ali.items()]
        #print(len(taxlist))
        #print(taxlist)
        datalib.phy2codeml(exp_dir + "multigene/" + "ppred" + basename + "_0_" + gene + ".ali.ali", simu_dir + gene + ".ali", original_ali)
        os.system("rm {0}ppred{1}_0_{2}.ali.ali".format(multi_dir, basename, gene))
        os.system("mv {0}ppred{1}_0_{2}.ali.trueparam {3}{2}.trueparam".format(multi_dir, basename, gene, simu_dir))
        os.system("mv {0}ppred{1}_0_{2}.ali.truesiteom {3}{2}.truesiteom".format(multi_dir, basename, gene, simu_dir))


