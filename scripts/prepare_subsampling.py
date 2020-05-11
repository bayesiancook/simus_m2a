#! /usr/bin/python3.5

import sys
import os
from runmm2a import * 
from random import sample
from ete3 import Tree

# input:
data_folder = sys.argv[1]
gene_file = sys.argv[2]
tree_file = sys.argv[3]
ngene = int(sys.argv[4])
nrep = int(sys.argv[5])
exp_folder = sys.argv[6]
basename = sys.argv[7]

data_dir = data_folder + "/"
exp_dir = exp_folder + "/"

machine = "p2chpd"
queue = "parallel"
path2bayescode = "/home_nfs/lbbe/nlartillot/bayescode/data/"
core = 16
nodes = 1
time = 24
mem = 16

if not os.path.exists(exp_folder):
    print("creating experiment folder")
    os.system("mkdir " + exp_folder)

os.system("cp " + tree_file + " " + exp_dir)

tree = Tree(tree_file)
taxlist = [leaf.name for leaf in tree]


with open(gene_file, 'r') as infile:
    gene_list = [gene.rstrip('\n').replace(".ali","") for gene in infile]

name2options = {"indmm2a" : "-x 1 1100 -nucrates ind -bl ind"}

included = dict()

for rep in range(nrep):

    list = sample(gene_list, ngene)
    for gene in list:
        included[gene] = 1

    with open(exp_dir + "sub{0}.list".format(rep), 'w') as outfile:
        outfile.write("{0}\n".format(ngene))
        for gene in list:
            outfile.write("{0}.ali\n".format(gene))

    for name in name2options:
        runname = "{0}{1}{2}".format(basename, name , rep)
        command = "multigenecodonm2a -d sub{0}.list -t {1} {2} {3}".format(rep, tree_file, name2options[name], runname)

        makebatch(command, runname, time=time, mem=mem, machine=machine, queue=queue, nodes=nodes, core=core, path2run=path2bayescode, path2batch = exp_dir)


list = [gene for gene in included]
make_single_gene_alignments(list, taxon_list = taxlist, all_taxa = True, from_dir = data_dir, to_dir = exp_dir)

