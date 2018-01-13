#! /usr/bin/env python

import sys, re, math
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

"""
def get_seq(reference, chr, start, end):

    seq = ""    
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        if item[0] == ">": continue
        seq = seq + item.rstrip('\n')

    if re.search(r'[^ACGTNacgtn]', seq) is not None:
        print >> sys.stderr, "The return value in get_seq function includes non-nucleotide characters:"
        print >> sys.stderr, seq
        sys.exit(1)

    return seq
"""

def get_seq(reference, chr, start, end):

    seq = ""
    for item in pysam.faidx(reference, chr + ":" + str(start) + "-" + str(end)):
        # if item[0] == ">": continue
        seq = seq + item.rstrip('\n')
    seq = seq.replace('>', '')
    seq = seq.replace(chr + ":" + str(start) + "-" + str(end), '')

    if re.search(r'[^ACGTNacgtn]', seq) is not None:
        print >> sys.stderr, "The return value in get_seq function includes non-nucleotide characters:"
        print >> sys.stderr, seq
        sys.exit(1)


    return seq


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


def generate_rss_pwm():

    # collected from Hesse et al., Gene. Dev., 1989
    heptamer_freq = [[0, 100, 0, 0], [99, 0, 1, 0], [0, 99, 0, 1], [87, 1, 3, 9], [10, 4, 82, 4], [3, 11, 1, 85], [13, 9, 76, 2]]
    nonamer_freq = [[72, 11, 9, 8], [5, 86, 3, 6], [83, 8, 6, 3], [73, 9, 8, 10], [91, 1, 3, 5], [97, 1, 0, 2], [87, 4, 2, 7], [4, 84, 4, 8], [5, 76, 4, 15]]
    background = [0.2, 0.3, 0.3, 0.2]

    heptamer_score = [[math.log(heptamer_freq[i][j] + 1) - math.log(sum(heptamer_freq[i]) + 4) - math.log(background[j]) for j in range(4)] for i in range(len(heptamer_freq))]
    nonamer_score = [[math.log(nonamer_freq[i][j] + 1) - math.log(sum(nonamer_freq[i]) + 4) - math.log(background[j]) for j in range(4)] for i in range(len(nonamer_freq))]

    return([heptamer_score, nonamer_score])


def get_pwm_score(seq, score_matrix):

    base2ind = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    score = 0
    for i in range(len(seq)):
        if seq[i] not in ['A', 'C', 'G', 'T']: continue
        score = score + score_matrix[i][base2ind[seq[i]]]

    return score


def get_max_rss_score(seq, pwm1, pwm2):

    pwm1_scores_plus, pwm2_scores_plus, pwm1_scores_minus, pwm2_scores_minus = [0] * len(seq), [0] * len(seq), [0] * len(seq), [0] * len(seq)

    rss_max_score = - float("inf")
    rss_max_pos = -1
    rss_max_spacer = -1
    rss_max_seq1 = ""
    rss_max_seq2 = ""
    rss_max_strand = "+"

    for i in range(len(seq)):
        if i + len(pwm1) > len(seq): continue
        pwm1_scores_plus[i] = get_pwm_score(seq[i:(i + len(pwm1))], pwm1) 

    for i in range(len(seq)):
        if i + len(pwm2) > len(seq): continue
        pwm2_scores_plus[i] = get_pwm_score(seq[i:(i + len(pwm2))], pwm2)

    for i in range(len(seq)):
        if i + len(pwm1) + len(pwm2) + 11 > len(seq): continue
        for j in [11, 12, 13, 22, 23, 24]:
            if i + len(pwm1) + len(pwm2) + j > len(seq): continue
            if pwm1_scores_plus[i] + pwm2_scores_plus[i + len(pwm1) + j] > rss_max_score:
                rss_max_score = pwm1_scores_plus[i] + pwm2_scores_plus[i + len(pwm1) + j]
                rss_max_pos = i
                rss_max_spacer = j
                rss_max_seq1 = seq[i:(i + len(pwm1))]
                rss_max_seq2 = seq[(i + len(pwm1) + j):(i + len(pwm1) + len(pwm2) + j)]
                rss_max_strand = '+'


    # for complementary strand
    seq_minus = reverse_complement(seq)
    for i in range(len(seq_minus)):
        if i + len(pwm1) > len(seq_minus): continue
        pwm1_scores_minus[i] = get_pwm_score(seq_minus[i:(i + len(pwm1))], pwm1)
 
    for i in range(len(seq_minus)):
        if i + len(pwm2) > len(seq_minus): continue
        pwm2_scores_minus[i] = get_pwm_score(seq_minus[i:(i + len(pwm2))], pwm2)

    for i in range(len(seq_minus)): 
        if i + len(pwm1) + len(pwm2) + 11 > len(seq_minus): continue
        for j in [11, 12, 13, 22, 23, 24]:
            if i + len(pwm1) + len(pwm2) + j > len(seq_minus): continue
            if pwm1_scores_minus[i] + pwm2_scores_minus[i + len(pwm1) + j] > rss_max_score:
                rss_max_score = pwm1_scores_minus[i] + pwm2_scores_minus[i + len(pwm1) + j]
                rss_max_pos = len(seq) - 1 - i
                rss_max_spacer = j
                rss_max_seq1 = seq_minus[i:(i + len(pwm1))]
                rss_max_seq2 = seq_minus[(i + len(pwm1) + j):(i + len(pwm1) + len(pwm2) + j)]
                rss_max_strand = '-'

    return [rss_max_score, rss_max_seq1, rss_max_seq2, rss_max_pos, rss_max_spacer, rss_max_strand]


