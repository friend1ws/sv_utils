#! /usr/bin/env python

import sys, subprocess
import my_seq
from header_info import *

def header_check(line):

    if "Chr_1" + '\t' + "Pos_1" + '\t' + "Dir_1" + '\t' + "Chr_2" + '\t' + "Pos_2" + '\t' + "Dir_2" + '\t' + "Inserted_Seq" + '\t' + "Variant_Type" in line:
        return True
    else:
        return False


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



def nonB_DB_dist_check(chr1, pos1, dir1, nonB_DB_tb, nonB_DB_type, margin = 500):

    nonB_DB_dist = "---"
 
    # check junction annotation for refGene 
    tabixErrorFlag = 0
    try:
        # this fetch procedure might be optimized more
        records = nonB_DB_tb.fetch(chr1, int(pos1) - margin * 2, int(pos1) + margin * 2)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            nonB_DB_dist_tmp = "---" 
            record = record_line.split('\t')
            if record[4] != nonB_DB_type: continue

            """
            motif_start = int(record[1]) + 1
            motif_end = int(record[2])
            if motif_start <= pos1 and pos1 <= motif_end:
                nonB_DB_dist_tmp = 0
            elif pos1 < motif_start and motif_start - pos1 <= margin:
                nonB_DB_dist_tmp = motif_start - pos1 if dir1 == "+" else - (motif_start - pos1)
            elif pos1 > motif_end and pos1 - motif_end <= margin:
                nonB_DB_dist_tmp = - (pos1 - motif_end) if dir1 == "+" else pos1 - motif_end
            """

            motif_center = (int(record[1]) + int(record[2]) + 1) / 2
            nonB_DB_dist_tmp = motif_center - pos1 if dir1 == "+" else pos1 - motif_center

            if abs(nonB_DB_dist_tmp) > margin: continue
            if nonB_DB_dist_tmp == "---": continue
            if nonB_DB_dist == "---" or abs(nonB_DB_dist_tmp) < abs(nonB_DB_dist): 
                nonB_DB_dist = nonB_DB_dist_tmp

    return nonB_DB_dist


