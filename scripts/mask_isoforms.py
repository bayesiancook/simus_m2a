#! /usr/bin/python3.5

import sys
import numpy as np
import os
import datalib
import translate

alphabet = "ACDEFGHIKLMNPQRSTVWY"
nstate = len(alphabet)

def mask_isoforms(aliname, outname, pseudocount = 0.05, cutoff = 0.45, width = 10):

    (ntaxa, nsite, nucali) = datalib.read_phylip(aliname)
    (ntaxa, naa, aaali) = translate.translate_nucali( (ntaxa, nsite, nucali) )

    # marginal entropy per site
    marg_ent = [0 for i in range(naa)]
    # entropy at each site conditional on state at previous site
    cond_ent = [0 for i in range(naa)]
    # corr = 1 - cond_ent / marg_ent
    corr = [0 for i in range(naa)]

    for i in range(1,naa):

        # initialize count matrix
        count_matrix = dict()
        for a in alphabet:
            row = dict()
            for b in alphabet:
                row[b] = pseudocount
            count_matrix[a] = row

        # fill-in count matrix
        for (tax,seq) in aaali.items():
            if (seq[i-1] in alphabet) and (seq[i] in alphabet):
                count_matrix[seq[i-1]][seq[i]] += 1

        # marginal entropy
        # marginal counts in second column regardless of first column
        marginal_count2 = [sum([count_matrix[a][b] for a in alphabet]) for b in alphabet]
        totcount = sum(marginal_count2)
        # renormalized as freqs
        marginal_freq2 = [count / totcount for count in marginal_count2]
        marg_ent[i] = -sum([p*np.log(p) for p in marginal_freq2])
    
        # conditional entropy
        marginal_count1 = [sum([count_matrix[a][b] for b in alphabet]) for a in alphabet]
        # matrix of conditional frequencies
        # row a: frequencies at 2d column, conditional on state at first column being a
        cond_freq_matrix = [[count_matrix[a][b] / marginal_count1[i] for b in alphabet] for i,a in enumerate(alphabet)]
        # row entropy
        row_entropy = [ -sum([p*np.log(p) for p in cond_freq_matrix[i]]) for i,a in enumerate(alphabet)]
        cond_ent[i] = sum([marginal_count1[i] / totcount * row_entropy[i] for i,a in enumerate(alphabet)])

        if marg_ent[i] > 1e-3:
            corr[i] = 1 - cond_ent[i] / marg_ent[i]

    mask = [int(max([corr[j] for j in range(max(0,i-width), min(naa,i+width+1))]) > cutoff) for i in range(naa)]
    excluded_nsite = 3*sum(mask)
    final_nsite = sum([1-m for m in mask])

    outali = dict()
    for (tax,seq) in nucali.items():
        outali[tax] = "".join([seq[3*i:3*(i+1)] for i in range(naa) if not mask[i]])

    with open(outname + ".ali", 'w') as outfile:
        outfile.write("{0}\t{1}\n".format(ntaxa,final_nsite))
        for (tax,seq) in outali.items():
            outfile.write("{0}  {1}\n".format(tax,seq))

    string_mask = "".join(["{}".format(m) for m in mask])
    with open(outname + ".mask", 'w') as outfile:
        outfile.write("{0}".format(string_mask))

    excluded_nsite = sum([1-m for m in mask])
    excluded_ali = dict()
    for (tax,seq) in aaali.items():
        excluded_ali[tax] = "".join([seq[i] for i in range(naa) if mask[i]])

    with open(outname + ".mask.fasta", 'w') as outfile:
        for (tax,seq) in excluded_ali.items():
            outfile.write(">{0}\n{1}\n".format(tax,seq))


aliname = sys.argv[1]
outname = sys.argv[2]

mask_isoforms(aliname, outname)

