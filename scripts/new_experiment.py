#! /usr/bin/python3.5

import sys
import os
from ete3 import Tree
from datalib import *

def new_experiment(data_folder, gene_file, taxon_file, tree_file, exp_folder):

    data_dir = data_folder + "/"
    exp_dir = exp_folder + "/"
    target_folder = exp_dir + "empirical"
    target_dir = target_folder + "/"

    if not os.path.exists(data_folder):
        print("data folder does not exist")
        raise

    if os.path.exists(exp_folder):
        print("experiment folder already exists")
        raise

    # get list of genes
    with open(gene_file, 'r') as file:
        gene_list = [gene.rstrip('\n').replace(".ali","") for gene in file]

    # get taxon list
    taxon_list = []
    if taxon_file != "alltaxa":
        with open(taxon_file) as file:
            taxon_list = [taxon.rstrip('\n') for taxon in file]
        ntaxa = len(taxon_list)
        print("number of taxa: {0}".format(ntaxa))

    # make target folder and data subfolder
    os.system("mkdir " + exp_folder)
    os.system("mkdir " + target_folder)
    os.system("mkdir " + target_dir + "data")

    ngene = len(gene_list)
    print("number of genes: {0}".format(ngene))

    # make single-gene alignments
    make_single_gene_alignments(gene_list, taxon_list = taxon_list, from_dir = data_dir, to_dir = target_dir + "data/")

    # make gene list file
    with open(target_dir + "all.list", 'w') as file:
        for gene in gene_list:
            file.write(gene + ".ali\n")

    # make tree file
    tree = Tree(tree_file)
    if taxon_list != []:
        tree.prune(taxon_list)
        tree.unroot()
    tree.write(outfile=target_dir + "all.tree",format=9)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("new_experiment.py data_folder gene_list taxon_list tree exp_folder")
        sys.exit()

    new_experiment(*sys.argv[1:])

