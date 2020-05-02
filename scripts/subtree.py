#! /beegfs/data/soft/anaconda3/bin/python3.7
from ete3 import Tree
import sys

if len(sys.argv) == 1:
    print("subtree.py treefile taxfile outfile")
    sys.exit()

treefilename = sys.argv[1]
taxfilename = sys.argv[2]
outfilename = sys.argv[3]

tree = Tree(treefilename)
with open(taxfilename, 'r') as taxfile:
    taxlist = [line.rstrip('\n') for line in taxfile]

tree.prune(taxlist)

tree.write(outfile=outfilename,format=9)



