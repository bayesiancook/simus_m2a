#! /usr/bin/python3.5
import sys
import os
from numpy import random

def genejack(list_name, ngene, nrep, outname):

    with open(list_name, 'r') as list_file:
        gene_list = [line.rstrip('\n') for line in list_file]

    tot_ngene = len(gene_list)
    if tot_ngene < nrep*ngene:
        print("error: gene list is too small")
        sys.exit()

    random.shuffle(gene_list)
    for rep in range(nrep):
        sub_list = gene_list[rep*ngene:(rep+1)*ngene]
        with open("{0}{1}.list".format(outname,rep), 'w') as outfile:
            for gene in sub_list:
                outfile.write(gene + "\n")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("genejack genelist ngene nrep outname")
        sys.exit()

    genejack(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])

