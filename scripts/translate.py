#! /usr/bin/python3.5

import sys
import numpy as np
import os
import datalib

univ_trans = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def translate_codon(codon):
    if codon in univ_trans:
        return univ_trans[codon]
    else:
        return "?"

def translate_nucseq(nucseq):
    nsite = len(nucseq)
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    return "".join([translate_codon(nucseq[3*i:3*(i+1)]) for i in range(naa)])

def translate_nucali(ali_triple):
    (ntaxa, nsite, nucali) = ali_triple
    if nsite % 3:
        print("error in translate: not multiple of 3")
        sys.exit()
    naa = nsite // 3
    aaali = dict()
    for (tax,nucseq) in nucali.items():
        aaali[tax] = translate_nucseq(nucseq)
    return (ntaxa, naa, aaali)


