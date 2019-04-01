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
        namelist_wocodeml = [name in namelist if name != "codeml"]
        method_gene_fdr(cutoff_list, namelist_wocodeml, score, trueposw, outname + "1")
        #gene_fdr(cutoff_list, namelist_wocodeml, score2, truepos2, outname + "2")
        #gene_fdr(cutoff_list, namelist_wocodeml, score3, truepos2, outname + "3")



