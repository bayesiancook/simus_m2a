
#! /usr/bin/python3.5

import sys
import numpy as np
import os
import datalib

def iso(aliname, outname, pseudocount = 0.05, cutoff = 0.45, width = 10):

    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    nstate = len(alphabet)

    # open alignment
    # fasta format with sequence on several lines
    ali = dict()
    with open(aliname, 'r') as alifile:

        taxname = ''
        ppseq = ''
        for line in alifile:
            line = line.rstrip('\n')
            if line[0] == '>':
                taxname = line[1:len(line)]
                if taxname != 'PP_CONS':
                    ali[taxname] = ''
            else:
                if taxname == 'PP_CONS':
                    ppseq = ppseq + line
                else:
                    ali[taxname] = ali[taxname] + line

        nsite = len(ppseq)

    # marginal entropy per site
    marg_ent = [0 for i in range(nsite)]
    # entropy at each site conditional on state at previous site
    cond_ent = [0 for i in range(nsite)]
    # corr = 1 - cond_ent / marg_ent
    corr = [0 for i in range(nsite)]

    for i in range(1,nsite):

        # initialize count matrix
        count_matrix = dict()
        for a in alphabet:
            row = dict()
            for b in alphabet:
                row[b] = pseudocount
            count_matrix[a] = row

        # fill-in count matrix
        for (tax,seq) in ali.items():
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

    # write alignment
    mask = [int(max([corr[j] for j in range(max(0,i-width), min(nsite,i+width+1))]) > cutoff) for i in range(nsite)]

    smooth = [np.mean([corr[j] for j in range(max(0,i-width), min(nsite,i+width+1))]) for i in range(nsite)]

    rec = dict()
    rec['N'] = 0
    rec['W'] = -0.5
    rec['R'] = -1.0
    recpp = [rec[p] for p in ppseq]

    # percentage of gaps
    gap_width = 10
    gap_score = [0 for i in range(nsite)]
    ent_score = [0 for i in range(nsite)]
    for i in range(nsite):
        jmin = max(0,i-gap_width)
        jmax = min(nsite,i+gap_width+1)

        ent_score[i] = np.mean([marg_ent[j] for j in range(jmin,jmax)])
        tot_npos = 0
        tot_ngap = 0
        for (tax,seq) in ali.items():
            seg = seq[jmin:jmax]
            ngap = sum([s not in alphabet for s in seg])
            if ngap != len(seg):
                tot_npos += len(seg)
                tot_ngap += ngap
        gap_score[i] = tot_ngap / tot_npos

    # stats of pos sel
    filtered_strong = sum([(mask[i] == 0) and (recpp[i] == -1.0) for i in range(nsite)])
    filtered_weak = sum([(mask[i] == 0) and (recpp[i] == -0.5) for i in range(nsite)])
    total_strong = sum([(recpp[i] == -1.0) for i in range(nsite)])
    total_weak = sum([(recpp[i] == -0.5) for i in range(nsite)])

    mean_gapscore_strong = 0
    if filtered_strong:
        mean_gapscore_strong = int(100 * np.mean([gap_score[i] for i in range(nsite) if (recpp[i] == -1.0) and (mask[i] == 0)]))
    mean_gapscore_weak  = 0
    if filtered_weak:
        mean_gapscore_weak = int(100 * np.mean([gap_score[i] for i in range(nsite) if (recpp[i] == -0.5) and (mask[i] == 0)]))
    mean_gapscore = int(100 * np.mean(gap_score))

    mean_entscore_strong = 0
    if filtered_strong:
        mean_entscore_strong = int(100 * np.mean([ent_score[i] for i in range(nsite) if (recpp[i] == -1.0) and (mask[i] == 0)]))
    mean_entscore_weak  = 0
    if filtered_weak:
        mean_entscore_weak = int(100 * np.mean([ent_score[i] for i in range(nsite) if (recpp[i] == -0.5) and (mask[i] == 0)]))
    mean_entscore = int(100 * np.mean(ent_score))

    print("raw", total_strong, total_weak)
    print("filtered", filtered_strong, filtered_weak)
    print("gap", mean_gapscore_strong, mean_gapscore_weak, mean_gapscore)
    print("entropy", mean_entscore_strong, mean_entscore_weak, mean_entscore)

    with open(outname + ".iso", 'w') as outfile:

        outfile.write("0\t0\t0\t0\0\0\n")
        for i in range(1,nsite):
            outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(corr[i], cond_ent[i], marg_ent[i], recpp[i], smooth[i], -0.9*mask[i], gap_score[i]))


aliname = sys.argv[1]
outname = sys.argv[2]
iso(aliname, outname)

