#! /usr/bin/python3.5

import sys
import os
from ete3 import Tree
from datalib import *
from makebatch import *

def prepare_subexperiment(exp_folder, prep_codeml = True, prep_single = True, prep_multi = True):

    exp_dir = exp_folder + "/"
    data_dir = exp_dir + "data/"

    # get gene list
    with open(exp_dir + "all.list") as file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in file]

    # get tree
    tree = Tree(exp_dir + "all.tree")

    # codeml

    if prep_codeml:
        # make target folder for codeml analyses
        # copy control files for codeml in this folder
        codeml_dir = exp_dir + "codeml/"
        if os.path.exists(codeml_dir):
            print("directory " + codeml_dir + " already exists")
            sys.exit()
        os.system("mkdir " + codeml_dir)

        # make single-gene alignments and trees, copy them into target folder
        make_single_gene_alignments(gene_list, from_dir = data_dir, to_dir = codeml_dir)
        make_single_gene_trees(gene_list, tree, from_dir = data_dir, to_dir = codeml_dir)

    # bayescode single-gene 

    if prep_single:

        # make directory
        single_dir = exp_dir + "singlegene/"
        if os.path.exists(single_dir):
            print("directory " + single_dir + " already exists")
            sys.exit()
        os.system("mkdir " + single_dir)

        # make single-gene alignments and trees, copy them into target folder
        make_single_gene_alignments(gene_list, from_dir = data_dir, to_dir = single_dir)
        make_single_gene_trees(gene_list, tree, from_dir = data_dir, to_dir = single_dir)

    # bayescode multi-gene 

    if prep_multi:
        # make directory
        multi_dir = exp_dir + "multigene/"
        if os.path.exists(multi_dir):
            print("directory " + multi_dir + " already exists")
            sys.exit()
        os.system("mkdir " + multi_dir)

        taxlist = [leaf.name for leaf in tree]
        # make single-gene alignments with all taxa, copy them into target folder
        make_single_gene_alignments(gene_list, taxon_list = taxlist, all_taxa = True, from_dir = data_dir, to_dir = multi_dir)

        # make list of genes + header in multigene folder
        with open(multi_dir + "all.list", 'w') as outlistfile:
            outlistfile.write("{0}\n".format(len(gene_list)))
            for gene in gene_list:
                outlistfile.write(gene + ".ali\n")

        # copy tree in multigene folder
        os.system("cp " + exp_dir + "all.tree " + multi_dir)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("prepare_subexperiment.py experiment/subexperiment")
        sys.exit()

    prepare_subexperiment(*sys.argv[1:])

