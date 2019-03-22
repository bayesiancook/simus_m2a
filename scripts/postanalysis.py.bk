#! /usr/bin/python3.5

import sys
import numpy
import os
import parsecodeml
import parsem2a
import parsemm2a
import simuparams
from scipy.stats import chi2

def gene_codeml_tdr(score, trueposw, outname):

    ngene = len(trueposw)
    print(ngene)

    with open(outname + ".genetdr_codeml", 'w') as outfile:
        etdr = 0
        tdr = 0
        n = 0
        for (gene,dlnl) in sorted(score.items(), key=lambda kv: kv[1], reverse=True):

            if dlnl>0:
                n = n + 1
                pval = (1 - chi2.cdf(2*dlnl,1))/2
                efdr = ngene*pval/n
                etdr = 1 - ngene*pval/n
                if trueposw[gene] > 0:
                    tdr = tdr + 1
                fdr = n - tdr
                    
                outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\n".format(gene,dlnl,pval,efdr,fdr/n))
                #outfile.write("{0:15s}\t{1:5.2f}\t{2:7.5f}\t{3:7.5f}\t{4:7.5f}\n".format(gene,dlnl,pval,etdr,tdr/n))

def gene_tdr(cutoff_list, namelist, score, trueposw, outname):

    methodgene_ndisc = dict()
    methodgene_tdr = dict()
    methodgene_etdr = dict()

    for name in namelist:
        if name != "codeml":
            gene_ndisc = dict()
            gene_tdr = dict()
            gene_etdr = dict()
            for cutoff in cutoff_list:
                gene_ndisc[cutoff] = 0
                gene_tdr[cutoff] = 0
                gene_etdr[cutoff] = 0

    for name in namelist:
        if name != "codeml":

            gene_ndisc = dict()
            gene_tdr = dict()
            gene_etdr = dict()
            for cutoff in cutoff_list:
                gene_ndisc[cutoff] = 0
                gene_tdr[cutoff] = 0
                gene_etdr[cutoff] = 0

            with open(outname + ".genetdr_" + name, 'w') as outfile:
                etdr = 0
                tdr = 0
                n = 0
                for (gene,pp) in sorted(score[name].items(), key=lambda kv: kv[1], reverse=True):
                    n = n + 1
                    etdr = etdr + pp
                    if trueposw[gene] > 0:
                        tdr = tdr + 1
                        
                    outfile.write("{0:15s}\t{1:5.2f}\t{2:5.2f}\t{3:5.2f}\n".format(gene,pp,etdr/n,tdr/n))

                    for cutoff in cutoff_list:
                        if pp > cutoff:
                            gene_ndisc[cutoff] = n
                            gene_tdr[cutoff] = tdr
                            gene_etdr[cutoff] = etdr

            methodgene_ndisc[name] = gene_ndisc
            methodgene_tdr[name] = gene_tdr
            methodgene_etdr[name] = gene_etdr

    with open(outname + ".genetdr", 'w') as outfile:
        outfile.write("")
        for name in namelist:
            if name != "codeml":
                outfile.write("  {0:>18s}".format(name))
        outfile.write("\n")

        outfile.write("cutoff")
        for name in namelist:
            if name != "codeml":
                outfile.write("  {0:>5s} {1:>5s} {2:>5s} ".format("disc", "etdr", "tdr"))
        outfile.write("\n")

        for cutoff in cutoff_list:

            outfile.write("{0:5.2f}".format(cutoff))
            for name in namelist:
                if name != "codeml":
                    tdr = 0
                    etdr = 0
                    ndisc = methodgene_ndisc[name][cutoff]
                    if ndisc:
                        tdr = methodgene_tdr[name][cutoff] / ndisc
                        etdr = methodgene_etdr[name][cutoff] / ndisc
                    outfile.write("  {0:5d} {1:5.2f} {2:5.2f} ".format(ndisc, etdr, tdr))

            outfile.write("\n")



def m2a_postanalysis(exp_folder, single_basename, multi_basename, outname = "m2a_postanalysis", burnin = 20, dlnlmin = 0, min_omega = 1.0):

    # compare codeml and m2a analyses

    # cutoff_list = [0.9]
    cutoff_list = [0.5, 0.7, 0.9]
    # cutoff_list = [0.5, 0.7, 0.9]

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
    selectedsites = dict()
    sitepp = dict()
    hyperparams = dict()

    namelist = []

    # parsing codeml results
    name = "codeml"
    codeml_res = parsecodeml.parse_list(codeml_dir, genelist)
    [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name]] = codeml_res[0:7]
    namelist.append(name)

    # parsing m2a results
    for name in single_basename:
        single_res = parsem2a.parse_list(name, genelist, int(burnin), path=single_dir, min_omega = min_omega)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name]] = single_res[0:9]
        namelist.append(name)

    # parsing mm2a results
    for name in multi_basename:
        multi_res = parsemm2a.parse_list(name, int(burnin), path=multi_dir)
        [score[name], posw[name], posom[name], minposom[name], maxposom[name], selectedsites[name], sitepp[name], score2[name], score3[name], hyperparams[name]] = multi_res[0:10]
        namelist.append(name)

    codeml_dlnl = score["codeml"]

    nsel = dict()
    etdr = dict()
    if fromsimu:
        tdr = dict()

    for cutoff in cutoff_list:
        cnsel = dict()
        cetdr = dict()
        if fromsimu:
            ctdr = dict()
        for name in namelist:
            gcnsel = dict()
            gcetdr = dict()
            if fromsimu:
                gctdr = dict()
            for gene in genelist:
                gcnsel[gene] = sum([(pp>cutoff) for (i,pp) in selectedsites[name][gene].items()])
                gcetdr[gene] = sum([pp for (i,pp) in selectedsites[name][gene].items() if pp > cutoff])
                if fromsimu:
                    gctdr[gene] = sum([ ((pp > cutoff) and (truesiteom[gene][i] > 1.0)) for (i,pp) in selectedsites[name][gene].items()])
            cnsel[name] = gcnsel
            cetdr[name] = gcetdr
            if fromsimu:
                ctdr[name] = gctdr
        nsel[cutoff] = cnsel
        etdr[cutoff] = cetdr
        if fromsimu:
            tdr[cutoff] = ctdr

    with open(outname + ".sortedparams", 'w') as geneoutfile:

        #
        # header 
        #

        geneoutfile.write("methods:\n")
        for i,name in enumerate(namelist):
            geneoutfile.write("{}\t{}\n".format(i,name))
        geneoutfile.write("\n")

        geneoutfile.write("{:15}".format("gene"))
        geneoutfile.write("  {:>6}".format("dlnl"))
        for i in range(len(namelist)-1):
            geneoutfile.write("  {:>4}{:1}".format("pp_",i+1))

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

            for cutoff in cutoff_list:
                geneoutfile.write(" ")
                for name in namelist:
                    geneoutfile.write("  {0:2d}".format(nsel[cutoff][name][gene]))
                    if fromsimu:
                        if name == "codeml":
                            geneoutfile.write(" ({0:2d})".format(tdr[cutoff][name][gene]))
                        else:
                            geneoutfile.write(" ({0:2d} ; {1:3.1f})".format(tdr[cutoff][name][gene], etdr[cutoff][name][gene]))

            geneoutfile.write("\n")

    with open(outname + ".sitepp", 'w') as siteoutfile:

        #
        # gene list in decreasing order of dlnl
        #

        for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):
            if lnl >= dlnlmin:
                for (i,pp) in selectedsites["codeml"][gene].items():
                    siteoutfile.write("{0}\t{1}\t{2}".format(gene,i+1,pp))
                    for name in namelist:
                        if name != "codeml":
                            siteoutfile.write("\t{0}".format(sitepp[name][gene][i]))
                    siteoutfile.write("\n")


    hypernamelist = ['purom_mean', 'purom_invconc', 'dposom_mean', 'dposom_invshape', 'purw_mean', 'purw_invconc', 'posw_mean', 'posw_invconc', 'pi']

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
        outfile.write("{0:10s}\t{1:^6s}\t{2:^6s}\t{3:^6s}\t{4:^6s}\t{5:^6s}\n".format("method", "pi", "posw", "stdev", "posom", "stedv"))
        if fromsimu:
            outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format("simu", pi, posw_mean, numpy.sqrt(posw_var), 1.0 + dposom_mean, numpy.sqrt(dposom_var)))

        #for name in namelist:
        #    (est_pi, est_purw_mean, est_purw_var, est_purw_invconc, est_posw_mean, est_posw_var, est_posw_invconc, est_purom_mean, est_purom_var, est_purom_invconc, est_posom_mean, est_posom_var, est_posom_invshape) = simuparams.get_empirical_hyperparams(posw[name], posw[name], posom[name], posom[name])
        #    outfile.write("{0:10s}\t{1:6.4f}\t{2:6.4f}\t{3:6.4f}\t{4:6.4f}\t{5:6.4f}\n".format(name, est_pi, est_posw_mean, numpy.sqrt(est_posw_var), est_posom_mean, numpy.sqrt(est_posom_var)))

    if fromsimu:

        gene_codeml_tdr(score["codeml"], trueposw, outname)

        cutoff_list = [0.5, 0.6, 0.7, 0.8, 0.9]

        methodsitendisc = dict()
        methodsitetdr = dict()
        methodsiteetdr = dict()

        for name in namelist:

            sitendisc = dict()
            sitetdr = dict()
            siteetdr = dict()

            for cutoff in cutoff_list:

                sitendisc[cutoff] = 0
                sitetdr[cutoff] = 0
                siteetdr[cutoff] = 0

                for (gene,lnl) in sorted(codeml_dlnl.items(), key=lambda kv: kv[1], reverse=True):
                    if lnl >= dlnlmin:

                        ndisc = sum([(pp > cutoff) for pp in sitepp[name][gene]])
                        etdr = sum([pp for pp in sitepp[name][gene] if pp > cutoff])
                        tdr = sum([ ((pp > cutoff) and (truesiteom[gene][i] > 1.0)) for i,pp in enumerate(sitepp[name][gene]) ])

                        sitendisc[cutoff] = sitendisc[cutoff] + ndisc
                        siteetdr[cutoff] = siteetdr[cutoff] + etdr
                        sitetdr[cutoff] = sitetdr[cutoff] + tdr

                    methodsitendisc[name] = sitendisc
                    methodsitetdr[name] = sitetdr
                    methodsiteetdr[name] = siteetdr
                

        with open(outname + ".sitetdr", 'w') as outfile:
            outfile.write("")
            for name in namelist:
                outfile.write("  {0:>18s}".format(name))
            outfile.write("\n")

            outfile.write("cutoff")
            for name in namelist:
                outfile.write("  {0:>5s} {1:>5s} {2:>5s} ".format("disc", "etdr", "tdr"))
            outfile.write("\n")

            for cutoff in cutoff_list:

                outfile.write("{0:5.2f}".format(cutoff))
                for name in namelist:

                    tdr = 0
                    etdr = 0
                    ndisc = methodsitendisc[name][cutoff]
                    if ndisc:
                        tdr = methodsitetdr[name][cutoff] / ndisc
                        etdr = methodsiteetdr[name][cutoff] / ndisc

                    outfile.write("  {0:5d} {1:5.2f} {2:5.2f} ".format(ndisc, etdr, tdr))

                outfile.write("\n")

        # gene true discovery rate

        cutoff_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        gene_tdr(cutoff_list, namelist, score, trueposw, outname + "1")
        gene_tdr(cutoff_list, namelist, score2, truepos2, outname + "2")
        gene_tdr(cutoff_list, namelist, score3, truepos2, outname + "3")
