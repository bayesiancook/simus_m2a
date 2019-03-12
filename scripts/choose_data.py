#! /usr/bin/python3.5
import sys
import os
import re
from ete3 import Tree

if len(sys.argv) == 1:
    print("choose_data.py data_folder list_name taxfile_name min_ntaxa min_nsite max_nsite target_folder")
    sys.exit()

data_folder = sys.argv[1]
listname = sys.argv[2]
taxfile_name = sys.argv[3]
min_ntaxa = int(sys.argv[4])
min_nsite = int(sys.argv[5])
max_nsite = int(sys.argv[6])
target_folder = sys.argv[7]

data_dir = data_folder + "/"
target_dir = target_folder + "/"

# make target directory
if os.path.exists(target_dir):
    print("directory " + target_dir + " already exists")
    sys.exit()
os.system("mkdir " + target_dir)

# get target taxon list
with open(taxfile_name, 'r') as taxfile:
    target_taxlist = [taxon.rstrip('\n') for taxon in taxfile]

# get gene list
with open(listname, 'r') as listfile:
    genelist = [gene.rstrip('\n').replace(".ali","") for gene in listfile]

all_missing = re.compile(r"^\?+$")

selection = []

for gene in genelist:

    # obtain subset of taxa from gene file
    with open(data_dir + gene + ".ali",'r') as genefile:
        (ntax,nsite) = genefile.readline().rstrip('\n').split()
        ali = dict()
        for line in genefile:
            entry = line.rstrip('\n').split()
            if not all_missing.match(entry[1]) and entry[0] in target_taxlist:
                ali[entry[0]] = entry[1]

        gene_taxlist = [taxon for taxon in ali]
        gene_ntaxa = len(gene_taxlist)
        gene_nsite = int(nsite)

        m = re.match(r".*_(ENSG\d+)_.*",gene) 
        if not m:
            print("error: gene does not match pattern")
            print(gene)
            raise
        gene_short_name = m.group(1)

        if gene_ntaxa >= min_ntaxa and gene_nsite >= min_nsite and gene_nsite <= max_nsite:
            # os.system("cp " + data_dir + gene + ".ali " + target_dir + "all" + gene_short_name + ".ali")
            selection.append(gene_short_name)
            with open(target_dir + gene_short_name + ".ali", 'w') as alifile:
                alifile.write("{0} {1}\n".format(gene_ntaxa, gene_nsite))
                for (taxon,seq) in ali.items():
                    alifile.write("{0}  {1}\n".format(taxon,seq))

with open(target_dir + "all.list", 'w') as listfile:
    for gene in selection:
        listfile.write(gene + ".ali\n")



