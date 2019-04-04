#! /usr/bin/python3.5
import sys
import os
import re
from ete3 import Tree

def multigene_choose_data(data_folder, listname, taxfile_name, target_folder, concat = True, min_ntaxa = 0, min_nsite = 0, max_nsite = 0):

    data_dir = data_folder + "/"
    target_dir = target_folder + "/"

    # make target directory
    if os.path.exists(target_dir):
        print("directory " + target_dir + " already exists")
        sys.exit()
    os.system("mkdir " + target_dir)

    # get target taxon list
    print("get taxon list", taxfile_name)
    with open(taxfile_name, 'r') as taxfile:
        target_taxlist = [taxon.rstrip('\n') for taxon in taxfile]

    print ("get gene list", listname)
    # get gene list
    with open(listname, 'r') as listfile:
        genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

    # all_missing = re.compile(r"^\?+$")
    # genus_name = re.compile(r"^([A-z][a-z]*)\_.*$")

    selection = dict()

    min = 0
    max = 0

    for gene in genelist:

        # obtain subset of taxa from gene file
        with open(data_dir + gene + ".ali",'r') as genefile:
            (ntax,nsite) = genefile.readline().rstrip('\n').split()
            ali = dict()
            for line in genefile:
                entry = line.rstrip('\n').split()

                # this is if we want to get only the genus name
                # m = re.match(r"^([A-Z][a-z]*)_.*$", entry[0])
                #if not m:
                #    print("error: taxon name does not match pattern")
                #    print(entry[0])
                #    sys.exit()
                #genus_name = m.group(1)
                #if not all_missing.match(entry[1]) and genus_name in target_taxlist:
                #    ali[genus_name] = entry[1]

                if entry[0] in target_taxlist:
                # to remove all missing taxa
                # if not all_missing.match(entry[1]) and entry[0] in target_taxlist:
                    ali[entry[0]] = entry[1]

            gene_taxlist = [taxon for (taxon,seq) in ali.items()]
            gene_ntaxa = len(gene_taxlist)
            gene_nsite = int(nsite)

            m = re.match(r"(ENSG\d{11})",gene) 
            if not m:
                print("error: gene does not match pattern")
                print(gene)
                sys.exit()

            gene_short_name = m.group(1)

            # print(gene_short_name, gene_ntaxa, gene_nsite)

            if gene_ntaxa >= min_ntaxa and gene_nsite >= min_nsite and (max_nsite == 0 or gene_nsite <= max_nsite):
                if min == 0 or min > gene_nsite :
                    min = gene_nsite
                if max < gene_nsite :
                    max = gene_nsite
                selection[gene_short_name] = (gene, gene_ntaxa, gene_nsite, ali)

    print("gene selection", len(selection))
    print("nsite min: ", min, " max: ", max)

    # output 

    # list of genes, with header
    with open(target_dir + "all.list", 'w') as listfile:

        listfile.write("{0}\n".format(len(selection)))

        for (gene_short_name, gene_full_name) in selection.items():
            listfile.write(gene_short_name + ".ali\n")

    if concat:
        with open(target_dir + "all.ali", 'w') as alifile:

            alifile.write("ALI\n")
            alifile.write("{0}\n".format(len(selection)))

            for (gene_short_name, gene_data) in selection.items():
                (gene_full_name, gene_ntaxa, gene_nsite, ali) = gene_data
                alifile.write("{0}\n".format(gene_short_name + ".ali"))
                alifile.write("{0} {1}\n".format(gene_ntaxa, gene_nsite))
                for taxon in target_taxlist:
                # for (taxon,seq) in ali.items():
                    if taxon in ali:
                        alifile.write("{0}  {1}\n".format(taxon,ali[taxon]))
                    else:
                        alifile.write("{0}  {1}\n".format(taxon,"?" * gene_nsite))

    else:

        for (gene_short_name, gene_data) in selection.items():
            with open(target_dir + gene_short_name + ".ali", 'w') as alifile:
                (gene_full_name, gene_ntaxa, gene_nsite, ali) = gene_data
                alifile.write("{0} {1}\n".format(gene_ntaxa, gene_nsite))
                for taxon in target_taxlist:
                # for (taxon,seq) in ali.items():
                    if taxon in ali:
                        alifile.write("{0}  {1}\n".format(taxon,ali[taxon]))
                    else:
                        alifile.write("{0}  {1}\n".format(taxon,"?" * gene_nsite))


if __name__ == "__main__":

    if len(sys.argv) == 1:
        print("multigene_choose_data.py data_folder list_name taxfile_name min_ntaxa min_nsite max_nsite target_folder")
        sys.exit()

    data_folder = sys.argv[1]
    listname = sys.argv[2]
    taxfile_name = sys.argv[3]
    min_ntaxa = int(sys.argv[4])
    min_nsite = int(sys.argv[5])
    max_nsite = int(sys.argv[6])
    target_folder = sys.argv[7]

    multigene_choose_data(data_folder, listname, taxfile_name, target_folder, concat = True, min_ntaxa = min_ntaxa, min_nsite = min_nsite, max_nsite = max_nsite)

