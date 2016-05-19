#! /usr/bin/env python

import pysam

def exact_alignment(seq1, seq2):

    A = [[0]* len(seq2) for i in range(len(seq1))]

    best_score = 0
    opt_coord = [0, 0]

    for i in range(len(seq1)):
        for j in range(len(seq2)):

            if i == 0:
                A[i][j] = 1 if seq1[i] == seq2[j] else 0
            elif j == 0:
                A[i][j] = 1 if seq1[i] == seq2[j] else 0
            else:
                A[i][j] = A[i - 1][j - 1] + 1 if A[i - 1][j - 1] > 0 and seq1[i] == seq2[j] else 0

            if A[i][j] > best_score:
                best_score = A[i][j]
                opt_coord = [i, j]

    opt_coord[0] = opt_coord[0] - best_score + 1
    opt_coord[1] = opt_coord[1] - best_score + 1
    return [best_score, opt_coord]

    """
    print '' + '\t' + '\t'.join(seq2)
    for i in range(len(seq1)):
        print seq1[i] + '\t' + '\t'.join([str(x) for x in A[i]])

    print ""
    print best_score
    print ""
    print opt_coord
    """


def get_seq(reference, chr, start, end):

    seq = ""    
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        if item[0] == ">": continue
        seq = seq + item.rstrip('\n')

    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))

