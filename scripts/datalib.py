import sys
import os

########

def read_phylip(infile_name):

    with open(infile_name, 'r') as infile:
        (ntax,npos) = infile.readline().rstrip('\n').split()
        ntaxa = int(ntax)
        nsite = int(npos)

        ali = dict()
        for line in infile:
            pair = line.rstrip('\n').split()
            if len(pair) != 2:
                print("error in readphylip: should have 2 fields per line")
                raise
            if len(pair[1]) != nsite:
                print("error in read_phylip: non matching number of sites")
                raise
            ali[pair[0]] = pair[1]
        if len(ali) != ntaxa:
            print("error in read_phylip: non matching number of taxa")
            raise

    return (ntaxa, nsite, ali)

########

# input : phylip alignment and (optional) taxon list
# output: codeml-compatible alignment (pruned from taxa not in taxon list if this list is provided)
# if all_taxa == True: add taxa (with all missing sequences) if not in original alignment

def phy2codeml(infile_name, outfile_name, taxon_list = [], all_taxa = False, force = False, prune_missing = False):

    if not force and os.path.exists(outfile_name):
        print("in phy2codeml: destination file already exists: ", outfile_name)
        raise

    (ntaxa, nsite, ali) = read_phylip(infile_name)

    with open(outfile_name, 'w') as outfile:

        all_missing  = "?" * nsite

        if all_taxa:
            if taxon_list == []:
                print("in phy2codeml: taxon list not provided and all_taxa == True")
                raise
            ntaxa2 = len(taxon_list)
            outfile.write("{0} {1}\n".format(ntaxa2,nsite))
            for tax in taxon_list:
                seq = all_missing
                if tax in ali:
                    seq = ali[tax]
                outfile.write("{0}  {1}\n".format(tax,seq))

        else:

            ntaxa2 = ntaxa;
            if taxon_list != []:
                ntaxa2 = len([tax for tax in ali if tax in taxon_list])

            outfile.write("{0} {1}\n".format(ntaxa2,nsite))
            for (tax,seq) in ali.items():
                if taxon_list == [] or tax in taxon_list:
                    if not prune_missing or seq != all_missing:
                        outfile.write("{0}  {1}\n".format(tax,seq))

########

# input : list of single gene alignments and taxon list
# output: list of single gene alignments compatible with codeml

def make_single_gene_alignments(gene_list, taxon_list = [], from_dir = "", to_dir = "", from_prefix = "", to_prefix = "", from_ext = ".ali", to_ext = ".ali", all_taxa = False, force = False, prune_missing  = False):

    if not os.path.exists(to_dir):
        os.system("mkdir " + to_dir)

    for gene in gene_list:
        phy2codeml(from_dir + from_prefix + gene + from_ext, to_dir + to_prefix + gene + to_ext, taxon_list = taxon_list, all_taxa = all_taxa, force = force, prune_missing  = prune_missing)


########

# input : tree and taxon list
# output: subtree spanned by taxa from taxon list

def make_subtree(tree, taxon_list):
    subtree = tree.copy()
    subtree.prune(taxon_list)
    subtree.unroot()
    return subtree

########

# input : list of single gene alignments and complete tree (and, optional, a list of taxa)
# output: list of single gene trees

def make_single_gene_trees(gene_list, tree, taxon_list = [], from_dir = "", to_dir = "", from_prefix = "", to_prefix = "", from_ext = ".ali", to_ext = ".tree", force = False):

    if not os.path.exists(to_dir):
        os.system("mkdir " + to_dir)

    for gene in gene_list:
        (ntaxa, nsite, ali) = read_phylip(from_dir + from_prefix + gene + from_ext)
        taxlist = [tax for tax in ali]
        if taxon_list != []:
            taxlist = [tax for tax in ali if tax in taxon_list]
        gene_tree = make_subtree(tree, taxlist)
        tree_outfile = to_dir + to_prefix + gene + to_ext
        if not force and os.path.exists(tree_outfile):
            print("in phy2codeml: destination file already exists: ", tree_outfile)
            raise
        gene_tree.write(outfile=tree_outfile,format=9)

########

