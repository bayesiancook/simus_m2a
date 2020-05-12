#! /usr/bin/python3.5

import sys
import numpy
import os

from fdr import * 
import parsecodeml
import parsem2a
import parsemm2a
import simuparams

# exp_folder: subexperiment
# single_basename: list of runs (under alternative priors) done with single-gene m2a
# multi_basename: list of runs (under alternative priors) done with multi-gene mm2a
# outname: base name for output
# <outname>.sortedparams  : main results tabulated across methods for all genes sorted by decreasing DlnL (codeml)
# <outname>.sitepp        : (deactivated) site BEB post probs
# <outname>.hyper         : median and 95CI for hyperparameters
# <outname>.genefdr       : posterior estimate of FDR

default_cutoff_list = [0.05, 0.1, 0.3]

def m2a_postanalysis(exp_folder, single_basename, multi_basename, outname = "m2a_postanalysis", single_burnin = 100, multi_burnin = 500, dlnlmin = 0, dlnlmax = 0, min_omega = 1.0, with_sites = False, refname = "indmm2a", genepp_cutoff = 0.5, sitepp_cutoff = 0.90, cutoff_list = default_cutoff_list):


    exp_dir = exp_folder + "/"
    codeml_dir = exp_dir + "codeml/"
    single_dir = exp_dir + "singlegene/"
    multi_dir = exp_dir + "multigene/"

    listname = "all.list"

    single_result = dict()

    # get gene list
    with open(exp_dir + listname, 'r') as listfile:
        genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

    ngene = len(genelist)

    # get number of taxa and sites per gene
    gene_ntaxa = dict()
    gene_nsite = dict()
    fromsimu = True
    for gene in genelist:
        with open(exp_dir + "data/" + gene + ".ali") as genefile:
            [ntaxa, nsite] = genefile.readline().rstrip('\n').split()
            gene_ntaxa[gene] = int(ntaxa)
            gene_nsite[gene] = int(nsite)
        if not os.path.exists(exp_dir + "data/" + gene + ".trueparam"):
            fromsimu = False
        if not os.path.exists(exp_dir + "data/" + gene + ".truesiteom"):
            fromsimu = False

    # get true param values (if from simu)
    # and calculate empirical means and variances from these true values
    if fromsimu:
        print("data are from simulation, getting true param values..")
        [truepurw, trueposw, truepurom, truedposom, truesiteom, truepos2] = simuparams.get_true_params(exp_folder)
        (pi, purw_mean, purw_var, purw_invconc, posw_mean, posw_var, posw_invconc, purom_mean, purom_var, purom_invconc, dposom_mean, dposom_var, dposom_invshape) = simuparams.get_empirical_hyperparams(truepurw, trueposw, truepurom, truedposom)

    truees = trueposw
    # truees = {gene : trueposw[gene]*dposom for (gene,dposom) in truedposom.items()}

    # for each method (codeml, m2a, mm2a, under any prior or settings)
    # and for all genes
    # we want posw and pp, posom and min max, and selected sites
    # posw[method_name][gene_name] ...

    posw = dict()
    score = dict()
    score2 = dict()
    score3 = dict()
    posom = dict()
    minposom = dict()
    maxposom = dict()
    meanom = dict()
    meanes = dict()
    selectedsites = dict()
    sitepp = dict()
    hyperparams = dict()

    namelist = ["codeml"] + single_basename + multi_basename;

    print("parsing codeml")
    # parsing codeml results
    name = "codeml"
    codeml_res = parsecodeml.parse_list(codeml_dir, genelist)
    [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name]] = codeml_res[0:7]
    meanes[name] = {gene : posw*posom[name][gene] for (gene,posw) in posw[name].items()}

    print("parsing single gene analyses")
    # parsing m2a results
    for name in single_basename:
        print(name)
        single_res = parsem2a.parse_list(name, genelist, single_burnin, path=single_dir, min_omega = min_omega)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name], meanom[name], meanes[name]] = single_res[0:11]

    print("parsing multi gene analyses")
    # parsing mm2a results
    for name in multi_basename:
        print(name)
        multi_res = parsemm2a.parse_list(name, multi_burnin, path=multi_dir, with_sites = with_sites)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name], hyperparams[name], meanes[name]] = multi_res[0:11]

    print("printing out sorted list in .sortedparams")

    codeml_dlnl = score["codeml"]

    with open(outname + ".sortedparams", 'w') as geneoutfile:

        # header 

        geneoutfile.write("methods:\n")
        for i,name in enumerate(namelist):
            geneoutfile.write("{}\t{}\n".format(i,name))
        geneoutfile.write("\n")

        geneoutfile.write("{:15}".format("gene"))
        geneoutfile.write("  {:>6}".format("dlnl"))
        for i in range(len(namelist)-1):
            geneoutfile.write("{:>7}{:1}".format("pp_",i+1))

        geneoutfile.write(" ")

        if fromsimu:
            geneoutfile.write("  {:>5}".format("w*"))

        for i in range(len(namelist)):
            geneoutfile.write("  {:>4}{:1}".format("w_",i))

        geneoutfile.write(" ")

        if fromsimu:
            geneoutfile.write("  {:>6}".format("om*"))

        for i in range(len(namelist)):
            geneoutfile.write("  {:>5}{:1}".format("om_",i))

        geneoutfile.write(" ")

        geneoutfile.write("\n")

        #
        # gene list in decreasing order of dlnl
        #

        for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):

            geneoutfile.write("{0:15s}".format(gene))
            for name in namelist:
                geneoutfile.write("  {0:6.2f}".format(score[name][gene]))

            geneoutfile.write(" ")

            if fromsimu:
                geneoutfile.write("  {0:5.3f}".format(trueposw[gene]))

            for name in namelist:
                geneoutfile.write("  {0:5.3f}".format(posw[name][gene]))

            geneoutfile.write(" ")

            if fromsimu:
                geneoutfile.write("  {0:6.2f}".format(1.0 + truedposom[gene]))

            for name in namelist:
                geneoutfile.write("  {0:6.2f}".format(posom[name][gene]))

            # site-level statistics -- deactivated

            #for cutoff in cutoff_list:
            #    geneoutfile.write(" ")
            #    for name in namelist:
            #        geneoutfile.write("  {0:2d}".format(nsel[cutoff][name][gene]))
            #        if fromsimu:
            #            if name == "codeml":
            #                geneoutfile.write(" ({0:2d})".format(tdr[cutoff][name][gene]))
            #            else:
            #                geneoutfile.write(" ({0:2d} ; {1:3.1f})".format(tdr[cutoff][name][gene], etdr[cutoff][name][gene]))
            geneoutfile.write("\n")

    #
    # site pps (against true status if available)
    #

    if with_sites:

        codeml_gene_set = []
        bayes_gene_set = []

        # for each gene
        # take all sites selected by codeml
        # sort their pp
        # for a given fdr cutoff: 0.95 / 0.90
        # count number of sites selected 
        # compute true fdr for each method
        # output: gene name, dlnl, p-value, gene pp, (#sites, #fdr) for codeml and bayes

        with open(outname + ".genesitepp", 'w') as outfile:
            for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):
                if dlnlmax == 0 or lnl <= dlnlmax:
                    if lnl >= dlnlmin:
                        codeml_gene_set.append(gene)
                    if score[refname][gene] > genepp_cutoff:
                        bayes_gene_set.append(gene)

            both_gene_set = [gene for gene in codeml_gene_set if gene in bayes_gene_set]
            codeml_not_bayes = [gene for gene in codeml_gene_set if gene not in bayes_gene_set]
            bayes_not_codeml = [gene for gene in bayes_gene_set if gene not in codeml_gene_set]

            outfile.write("codeml cutoff: {0}\n".format(dlnlmin))
            outfile.write("bayes cutoff: {0}\n".format(genepp_cutoff))
            outfile.write("#genes by codeml: {0}\n".format(len(codeml_gene_set)))
            outfile.write("#genes by bayes: {0}\n".format(len(bayes_gene_set)))
            outfile.write("#genes by both: {0}\n".format(len(both_gene_set)))
            outfile.write("\n")

            # number of sites selected by codeml
            # fdr over those sites
            # mean codeml site pp over those sites
            c_n = 0
            c_fdr = 0
            c_meansitepp = 0

            # number of sites selected by bayes
            # fdr over those sites
            # mean bayes site pp over those sites
            b_n = 0
            b_fdr = 0
            b_meansitepp = 0

            # number of sites selected by codeml but not by bayes
            # fdr over those sites
            # mean codeml site pp over those sites
            # mean bayes site pp over those sites
            cb_n = 0
            cb_fdr = 0
            cbc_meansitepp = 0
            cbb_meansitepp = 0

            # number of sites selected by bayes but not by codeml
            # fdr over those sites
            # mean codeml site pp over those sites
            # mean bayes site pp over those sites
            bc_n = 0
            bc_fdr = 0
            bcc_meansitepp = 0
            bcb_meansitepp = 0

            for gene in both_gene_set:
                for (i,c_pp) in selectedsites["codeml"][gene].items():

                    b_pp = sitepp[refname][gene][i]
                    notpos = truesiteom[gene][i] == 1.0

                    if c_pp > sitepp_cutoff:
                        c_n = c_n + 1
                        c_meansitepp = c_meansitepp + c_pp
                        if notpos:
                            c_fdr = c_fdr + 1

                    if b_pp > sitepp_cutoff:
                        b_n = b_n + 1
                        b_meansitepp = b_meansitepp + b_pp
                        if notpos:
                            b_fdr = b_fdr + 1

                    if c_pp > sitepp_cutoff and b_pp <= sitepp_cutoff:
                        cb_n = cb_n + 1
                        cbc_meansitepp = cbc_meansitepp + c_pp
                        cbb_meansitepp = cbb_meansitepp + b_pp
                        if notpos:
                            cb_fdr = cb_fdr + 1
                        
                    if b_pp > sitepp_cutoff and c_pp <= sitepp_cutoff:
                        bc_n = bc_n + 1
                        bcc_meansitepp = bcc_meansitepp + c_pp
                        bcb_meansitepp = bcb_meansitepp + b_pp
                        if notpos:
                            bc_fdr = bc_fdr + 1

            c_meansitepp = c_meansitepp / c_n
            c_fdr = c_fdr / c_n

            b_meansitepp = b_meansitepp / b_n
            b_fdr = b_fdr / b_n

            if cb_n:
                cbc_meansitepp = cbc_meansitepp / cb_n
                cbb_meansitepp = cbb_meansitepp / cb_n
                cb_fdr = cb_fdr / cb_n

            if bc_n:
                bcc_meansitepp = bcc_meansitepp / bc_n
                bcb_meansitepp = bcb_meansitepp / bc_n
                bc_fdr = bc_fdr / bc_n

            outfile.write("{0} sites found by codeml\n".format(c_n))
            outfile.write("fdr : {0}\n".format(c_fdr))
            outfile.write("codeml e-fdr : {0}\n".format(1-c_meansitepp))
            outfile.write("\n")

            outfile.write("{0} sites found by bayes\n".format(b_n))
            outfile.write("fdr : {0}\n".format(b_fdr))
            outfile.write("bayes e-fdr : {0}\n".format(1-b_meansitepp))
            outfile.write("\n")

            outfile.write("{0} sites found by codeml but not bayes\n".format(cb_n))
            outfile.write("fdr : {0}\n".format(cb_fdr))
            outfile.write("codeml e-fdr : {0}\n".format(1-cbc_meansitepp))
            outfile.write("bayes e-fdr : {0}\n".format(1-cbb_meansitepp))
            outfile.write("\n")

            outfile.write("{0} sites found by bayes but not codeml\n".format(bc_n))
            outfile.write("fdr : {0}\n".format(bc_fdr))
            outfile.write("bayes e-fdr : {0}\n".format(1-bcb_meansitepp))
            outfile.write("codeml e-fdr : {0}\n".format(1-bcc_meansitepp))
            outfile.write("\n")

        with open(outname + ".compsitepp", 'w') as siteoutfile:
            siteoutfile.write("{0:15s}{1:5s}".format("gene", "site"))
            for name in namelist:
                siteoutfile.write("    {0:6s}".format("g_" + name))
            if fromsimu:
                siteoutfile.write("   {0:5s}".format("trueom"))
            for name in namelist:
                siteoutfile.write("   {0:5s}".format("s_" + name))
            siteoutfile.write("\n")

            for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):
                if score[refname][gene] > genepp_cutoff:
                # if lnl >= dlnlmin:

                    for (i,pp) in selectedsites["codeml"][gene].items():
                        siteoutfile.write("{0:15s}".format(gene))
                        siteoutfile.write("{0:5d}".format(i+1))
                        for name in namelist:
                            siteoutfile.write("   {0:6.2f}".format(score[name][gene]))
                        if fromsimu:
                            siteoutfile.write("  {0:5.2f}".format(truesiteom[gene][i]))
                        siteoutfile.write("  {0:5.2f}".format(pp))
                        for name in namelist:
                            if name != "codeml":
                                siteoutfile.write("  {0:5.2f}".format(sitepp[name][gene][i]))
                        siteoutfile.write("\n")


    hypernamelist = ['purom_mean', 'purom_invconc', 'dposom_mean', 'dposom_invshape', 'purw_mean', 'purw_invconc', 'posw_mean', 'posw_invconc', 'pi']

    print("printing out estimated hyperparameters")
    with open(outname + ".hyper", 'w') as outfile:

        # posterior mean and credible interval for hyperparams of multigene runs
        outfile.write("{0:20s}".format("param"))
        for name in hyperparams:
            outfile.write("\t{0:^24s}".format(name))
        outfile.write('\n')
        for i,name in enumerate(hypernamelist):
            outfile.write("{0:20s}".format(name))
            for name in hyperparams:
                (mean,min,max) = hyperparams[name][i]
                outfile.write("\t{0:6.4f} ({1:6.4f} , {2:6.4f})".format(mean,min,max))
            outfile.write('\n')

        outfile.write('\n')

        # mean and stdev of distribution of true versus post mean estimates for the 4 key parameters
        # outfile.write("{0:10s}\t{1:^6s}\t{2:^6s}\t{3:^6s}\t{4:^6s}\t{5:^6s}\n".format("method", "pi", "posw", "stdev", "posom", "stedv"))
        # if fromsimu:
        #    outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format("simu", pi, posw_mean, numpy.sqrt(posw_var), 1.0 + dposom_mean, numpy.sqrt(dposom_var)))

        #for name in namelist:
        #    (est_pi, est_purw_mean, est_purw_var, est_purw_invconc, est_posw_mean, est_posw_var, est_posw_invconc, est_purom_mean, est_purom_var, est_purom_invconc, est_posom_mean, est_posom_var, est_posom_invshape) = simuparams.get_empirical_hyperparams(posw[name], posw[name], posom[name], posom[name])
        #    outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format(name, est_pi, est_posw_mean, numpy.sqrt(est_posw_var), est_posom_mean, numpy.sqrt(est_posom_var)))

    print("gene-level fdr")
    truepos = dict()
    if fromsimu:
        truepos = trueposw

    return method_gene_fdr(cutoff_list, namelist, score, posw, trueposw, meanes, truees, gene_nsite, outname)


def full_m2a_postanalysis(exp_folder, simu_list, single_basename, multi_basename, outname, single_burnin = 100, multi_burnin = 500, with_sites = False, fdr_cutoff_list = default_cutoff_list, fields = ["ndisc", "fdr", "e-fnr", "fnr"]):

    simu_ret = dict()

    exp_dir = exp_folder + "/"
    res_dir = exp_dir + "results/"

    for simu in simu_list:
        print(simu)
        simu_ret[simu] = m2a_postanalysis(exp_dir + simu, single_basename, multi_basename, single_burnin = single_burnin, multi_burnin = multi_burnin, with_sites = with_sites, cutoff_list = fdr_cutoff_list, outname = res_dir + simu)

    nf = len(fields)
    namelist = ["df2_codeml"] + single_basename + multi_basename
    # namelist = ["df1_codeml", "df2_codeml", "mixdf1_codeml"] + single_basename + multi_basename

    with open(res_dir + outname + ".summary", 'w') as outfile:

        outfile.write("{0:>18s}".format(""))
        for cutoff in fdr_cutoff_list:
            outfile.write('{1:^{0}.2f}'.format(6*nf,cutoff))
        outfile.write("\n")
        outfile.write("\n")

        for simu in simu_list:

            outfile.write("{0:>18s}".format(simu))
            for cutoff in fdr_cutoff_list:
                for field in fields:
                    outfile.write(" {0:>5s}".format(field))
            outfile.write("\n")

            for name in namelist:

                outfile.write("{0:>18s}".format(name))

                for cutoff in fdr_cutoff_list:
                    for field in fields:
                        if (field == "n"):
                            outfile.write(" {0:5d}".format(simu_ret[simu][field][name][cutoff]))
                        else:
                            if name in simu_ret[simu][field]:
                                #if field in ["es", "ees"]:
                                #    outfile.write(" {0:5.0f}".format(simu_ret[simu][field][name][cutoff]))
                                #else:
                                #    outfile.write(" {0:5.2f}".format(simu_ret[simu][field][name][cutoff]))
                                outfile.write(" {0:5.2f}".format(simu_ret[simu][field][name][cutoff]))
                            else:
                                outfile.write(" {0:^5s}".format("-"))
                        
                outfile.write("\n")

            outfile.write("\n")

