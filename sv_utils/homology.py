#! /usr/bin/env python

import my_seq


def check_homology(chr1, pos1, dir1, chr2, pos2, dir2, reference, seq_size, min_match_size = 2, match_margin_ratio_thres = 0.5):

    seq1_out = my_seq.get_seq(reference, chr1, int(pos1) + 1, int(pos1) + seq_size)
    seq1_in = my_seq.reverse_complement(my_seq.get_seq(reference, chr1, int(pos1) - seq_size + 1, int(pos1)))
    if dir1 == "-": seq1_out, seq1_in = seq1_in, seq1_out

    seq2_out = my_seq.get_seq(reference, chr2, int(pos2) + 1, int(pos2) + seq_size)
    seq2_in = my_seq.reverse_complement(my_seq.get_seq(reference, chr2, int(pos2) - seq_size + 1, int(pos2)))
    if dir2 == "-": seq2_out, seq2_in = seq2_in, seq2_out

    match1, coord1 = my_seq.exact_alignment(seq1_out, seq2_in)
    match2, coord2 = my_seq.exact_alignment(seq2_out, seq1_in)

    if max(coord1) > match_margin_ratio_thres * match1: match1 = 0
    if max(coord2) > match_margin_ratio_thres * match2: match2 = 0
    
    match = max(match1, match2) if max(match1, match2) >= min_match_size else 0
    return match


