#! /usr/bin/python3.5

import sys
import numpy
import os


def make_table(inname, outname, simu_list = [], trans_simu = dict(), method_list = [], trans_method = dict(), entry_list = [], trans_entry = dict(), target = "target FDR"):

    # read infile
    with open(inname, 'r') as infile:
        cutoff_list = infile.readline().rstrip('\n').split()
        print(cutoff_list)
        
        newsimu = False
        nentries = 0
        entries = []

        content = dict()
        methods = []
        simus = []

        for line in infile:
            line.rstrip('\n')
            fields = line.split()
            if not fields:
                newsimu = True
            else:
                if newsimu:
                    simuname = fields[0]
                    simus.append(simuname)
                    content[simuname] = dict()
                    for cutoff in cutoff_list:
                        content[simuname][cutoff] = dict()

                    k = (len(fields)-1) // len(cutoff_list)
                    if k*len(cutoff_list) != len(fields)-1:
                        print("parsing error: number of fields not multiple of number of cutoff thresholds")
                        sys.exit()
                    if not nentries:
                        nentries = k
                    else:
                        if nentries != k:
                            print("parsing error: incorrect number of fields")
                            sys.exit()

                    if not entries:
                       entries = fields[1:k+1]

                    for cutoff in cutoff_list:
                       for entry in entries:
                           content[simuname][cutoff][entry] = dict()

                    newsimu = False

                else:
                    if len(fields) != 1 + nentries * len(cutoff_list):
                        print("error when reading values: incorrect number of fields\n")
                        sys.exit()

                    method = fields[0]
                    if len(simus) == 1:
                        methods.append(method)
                    else:
                        if method not in methods:
                            print("error: did not find method")
                            print("{0} for simu {1}".format(method, simuname))
                            sys.exit()
                    c = 0
                    k = 0
                    for val in fields[1:]:
                        content[simuname][cutoff_list[c]][entries[k]][method] = val
                        k = k+1
                        if k == nentries:
                            c = c+1
                            k = 0

    if not entry_list:
        entry_list = entries
    else:
        for entry in entry_list:
            if entry not in entries:
                print("error: did not find entry {0}".format(entry))
                sys.exit()
    if not trans_entry:
        trans_entry = {entry:entry for entry in entry_list}
    else:
        for entry in entry_list:
            if entry not in trans_entry:
                print("error: did not find entry {0} in translation table".format(entry))
                sys.exit()

    if not method_list:
        method_list = methods
    else:
        for method in method_list:
            if method not in methods:
                print("error: did not find method {0}".format(method))
                sys.exit()
    if not trans_method:
        trans_method = {method:method for method in method_list}
    else:
        for method in method_list:
            if method not in trans_method:
                print("error: did not find method {0} in translation table".format(method))
                sys.exit()

    if not simu_list:
        simu_list = simus
    else:
        for simu in simu_list:
            if simu not in simus:
                print("error: did not find simu {0}".format(simu))
                sys.exit()
    if not trans_simu:
        trans_simu = {simu:simu for simu in simu_list}
    else:
        for simu in simu_list:
            if simu not in trans_simu:
                print("error: did not find simu {0} in translation table".format(simu))
                sys.exit()

    subnentries = len(trans_entry)

    with open(outname + ".tex", 'w') as outfile:

        tabul = 'r' + 'c' * subnentries * len(cutoff_list)
        outfile.write("\\begin{{tabular}}{{{0}}}\n".format(tabul))

        outfile.write(" & \\multicolumn{{ {0} }}{{c}}{{ {1} }}".format(subnentries*len(cutoff_list), target))
        outfile.write(r'\\')
        outfile.write("\n")

        for cutoff in cutoff_list:
            outfile.write("& \\multicolumn{{ {0} }}{{c}}{{${1}$}}".format(subnentries,cutoff))
        outfile.write(r'\\')
        outfile.write("\n")

        for simu in simu_list:
            simutex = trans_simu[simu]
            outfile.write("{0}".format(simutex))
            for cutoff in cutoff_list:
                for entry in entry_list:
                    entrytex = trans_entry[entry]
                    outfile.write("& {0}".format(entrytex))
            outfile.write(r'\\')
            outfile.write("\n")

            # outfile.write("\\hline\n")

            for method in method_list:
                methodtex = trans_method[method]
                outfile.write("{0:>18s}".format(trans_method[method]))

                for cutoff in cutoff_list:
                    for entry in entry_list:
                        entrytex = trans_entry[entry]
                        outfile.write("& ${0}$".format(content[simu][cutoff][entry][method]))
                outfile.write(r'\\')
                outfile.write("\n")

            outfile.write("\\hline\n")
            # outfile.write(r'\\')
            # outfile.write("\n")

        outfile.write(r'\end{tabular}{}')
        outfile.write("\n")



if __name__ == "__main__":


    print("table 2")

    target = "target FDR"

    inname = "allsimu.summary"
    outname = "table2"

    trans_simu = {
            "simu30" : "\\bf sim 30\\%",
            "simu10" : "\\bf sim 10\\%",
            "simu03" : "\\bf sim 03\\%",
            "simu30_shrink_dposom03" : "\\bf sim 30\\% smaller $\\omega_+$",
            "simu30_shrink_posw03" : "\\bf sim 30\\% smaller $w_+$"
            }

    trans_method = {
            "df2_codeml" : "codeml",
            "unconsindmm2a" : "bayes",
            "indmm2a" : "bayes modal",
            "unconsshrunkenmm2a" : "bayes shrunken",
            "shrunkenmm2a" : "bayes shrunken modal"
            }

    trans_entry = {
            "n" : "n",
            "fdr" : "fdr",
            "efnr" : "efnr",
            "fnr" : "fnr"
            }

    simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_posw03", "simu30_shrink_dposom03"]
    method_list = ["df2_codeml", "unconsindmm2a", "indmm2a"]
    entry_list = ["n", "fdr", "efnr", "fnr"]

    make_table(inname, outname, simu_list = simu_list, trans_simu = trans_simu, method_list = method_list, trans_method = trans_method, entry_list = entry_list, trans_entry = trans_entry, target = target)

    print("table 3")

    target = "target FDR"

    inname = "allsimu.summary"
    outname = "table3"

    trans_simu = {
            "simu10" : "\\bf sim independent",
            "simushrunken" : "\\bf sim shrunken",
            "simushared" : "\\bf sim shared"
            }

    trans_method = {
            "indmm2a" : "bayes independent",
            "shrunkenmm2a" : "bayes shrunken",
            "sharedmm2a" : "bayes shared"
            }


    simu_list = ["simu10", "simushrunken", "simushared"]
    method_list = ["indmm2a", "shrunkenmm2a", "sharedmm2a"]
    entry_list = ["n", "fdr", "efnr", "fnr"]

    make_table(inname, outname, simu_list = simu_list, trans_simu = trans_simu, method_list = method_list, trans_method = trans_method, entry_list = entry_list, trans_entry = trans_entry, target = target)

    print("table 4")

    inname = "allsimupos.es_summary"
    outname = "table4"

    target = "target \% sites under pos. sel."

    trans_simu = {
            "simu30" : "\\bf sim 30\\%",
            "simu10" : "\\bf sim 10\\%",
            "simu03" : "\\bf sim 03\\%",
            "simu30_shrink_dposom03" : "\\bf sim 30\\% smaller $\\omega_+$",
            "simu30_shrink_posw03" : "\\bf sim 30\\% smaller $w_+$",
            "simushrunken" : "\\bf sim shrunken"
            }

    trans_method = {
            "df2_codeml" : "codeml",
            "indmm2a" : "bayes modal",
            "unconsindmm2a" : "bayes"
            }

    trans_entry = {
            "n" : "n",
            "%es" : "\%sites",
            "efdr" : "e-fdr",
            "fdr" : "fdr"
            }

    simu_list = ["simushrunken", "simu30", "simu10", "simu03", "simu30_shrink_posw03", "simu30_shrink_dposom03"]
    method_list = ["unconsindmm2a", "indmm2a"]
    # method_list = ["df2_codeml", "unconsindmm2a", "indmm2a"]
    entry_list = ["n", "%es", "efdr", "fdr"]

    make_table(inname, outname, simu_list = simu_list, trans_simu = trans_simu, method_list = method_list, trans_method = trans_method, entry_list = entry_list, trans_entry = trans_entry, target = target)


    print("table 5")

    target = "target \% excess $dN/dS$ ($\\Delta \\omega_+$)"

    inname = "allsimuexcess.es_summary"
    outname = "table5"

    trans_simu = {
            "simu30" : "\\bf sim 30\\%",
            "simu10" : "\\bf sim 10\\%",
            "simu03" : "\\bf sim 03\\%",
            "simu30_shrink_dposom03" : "\\bf sim 30\\% smaller $\\omega_+$",
            "simu30_shrink_posw03" : "\\bf sim 30\\% smaller $w_+$"
            }

    trans_method = {
            "indmm2a" : "bayes modal",
            "unconsindmm2a" : "bayes"
            }

    trans_entry = {
            "n" : "n",
            "%es" : "\%$\\Delta \\omega_+$",
            "efdr" : "e-fdr",
            "fdr" : "fdr"
            }

    simu_list = ["simu30", "simu10", "simu03", "simu30_shrink_posw03", "simu30_shrink_dposom03"]
    method_list = ["unconsindmm2a", "indmm2a"]
    entry_list = ["n", "%es", "efdr", "fdr"]

    make_table(inname, outname, simu_list = simu_list, trans_simu = trans_simu, method_list = method_list, trans_method = trans_method, entry_list = entry_list, trans_entry = trans_entry, target = target)

