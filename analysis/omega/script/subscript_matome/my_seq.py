#! /usr/bin/env python

import pysam

def get_seq(reference, pos):

    seq = ""    
    for item in pysam.faidx(reference, pos):
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(pos, '')

    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


