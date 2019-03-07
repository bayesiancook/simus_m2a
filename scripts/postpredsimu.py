import sys
import os
import re

def single_gene_postpred(exp_folder, listname, basename, target_folder, lnl_min = -1.0, lnl_max = 0.0, burnin = 10, every = 10, shrink_posw = 1.0, shrink_dposom = 1.0, nrep = 1):

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

    # make ppred simulations for all genes
    for (gene,lnl) in sorted(genelist.items(), key=lambda kv: kv[1], reverse=True):

        # check trees are there
        if not os.path.exists(single_dir + gene + ".tree"):
            print("gene tree " + single_dir + gene + ".tree does not exist")
            raise
        os.system("cp " + single_dir + gene + ".tree " + simu_dir)

        # readcodonm2a with ppred option
        ppred_command = "readcodonm2a -x {0} {1} {2} -ppred -shrinkposw {3} -shrinkdposom {4} {5}{6}".format(burnin, every, until, shrink_posw, shrink_dposom, basename, gene)
        if ((lnl_max > 0) and (lnl > lnl_max)) or (lnl < lnl_min):
            ppred_command = "readcodonm2a -x {0} {1} {2} -ppred -null {3}{4}".format(burnin,every,until,basename,gene)

        # execute ppred and move resulting files in simu folder
        os.chdir(single_dir)
        os.system(ppred_command)
        os.chdir("../../")
        os.system("mv {0}ppred{1}{2}_{3}.ali {4}{2}.ali".format(single_dir, basename, gene, rep, simu_dir))
        os.system("mv {0}ppred{1}{2}_{3}.trueparam {4}{2}.trueparam".format(single_dir, basename, gene, rep, simu_dir))
        os.system("mv {0}ppred{1}{2}_{3}.truesiteom {4}{2}.truesiteom".format(single_dir, basename, gene, rep, simu_dir))

