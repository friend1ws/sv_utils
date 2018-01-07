#! /usr/bin/env python

import sys, subprocess
import my_seq
from header_info import *

def header_check(line):

    if "Chr_1" + '\t' + "Pos_1" + '\t' + "Dir_1" + '\t' + "Chr_2" + '\t' + "Pos_2" + '\t' + "Dir_2" + '\t' + "Inserted_Seq" + '\t' + "Variant_Type" in line:
        return True
    else:
        return False


# get the closest exon and distances to them
def distance_to_closest(chr1, pos1, chr2, pos2, ref_tb, is_exon = True, search_max = 500000):

    cur_dist = search_max
    target_gene = []

    # check junction annotation for refGene for the first break point
    tabixErrorFlag = 0
    try:
        search_start = max(0, int(pos1) - search_max)
        records = ref_tb.fetch(chr1, search_start, int(pos1) + search_max)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if is_exon == False and record[4] != "coding": continue 
            # if within exon
            temp_dist = search_max
            if int(record[1]) < int(pos1) and int(pos1) <= int(record[2]):
                temp_dist = 0
            else:
                temp_dist = max(int(record[1]) - int(pos1) + 1, int(pos1) - int(record[2]))

            if temp_dist < cur_dist:
                target_gene = [record[3]]
                cur_dist = temp_dist
            elif temp_dist == cur_dist:
                target_gene.append(record[3])


    # check junction annotation for refGene for the second break point
    tabixErrorFlag = 0
    try:
        search_start = max(0, int(pos2) - search_max)
        records = ref_tb.fetch(chr2, search_start, int(pos2) + search_max)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
    
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if is_exon == False and record[4] != "coding": continue
            # if within exon
            temp_dist = search_max
            if int(record[1]) < int(pos2) and int(pos2) <= int(record[2]):
                temp_dist = 0
            else:
                temp_dist = max(int(record[1]) - int(pos2) + 1, int(pos2) - int(record[2]))
                
            if temp_dist < cur_dist:
                target_gene = [record[3]]
                cur_dist = temp_dist
            elif temp_dist == cur_dist:
                target_gene.append(record[3])


    target_gene = list(set(target_gene))
    return [cur_dist, target_gene]


def get_gene_annotation(chr1, pos1, chr2, pos2, gene_tb, exon_tb):

    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = gene_tb.fetch(chr1, int(pos1) - 1, int(pos1))
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    gene1 = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            gene1.append(record[3])

    if len(gene1) == 0: gene1.append("---")
    gene1 = list(set(gene1))
    ##########

    ##########
    # check gene annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = gene_tb.fetch(chr2, int(pos2) - 1, int(pos2))
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    gene2 = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            gene2.append(record[3])

    if len(gene2) == 0: gene2.append("---")
    gene2 = list(set(gene2))
    ##########

    ##########
    # check exon annotation for the side 1
    tabixErrorFlag = 0
    try:
        records = exon_tb.fetch(chr1, int(pos1) - 1, int(pos1))
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    exon1 = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            exon1.append(record[3])

    if len(exon1) == 0: exon1.append("---")
    exon1 = list(set(exon1))
    ##########

    ##########
    # check exon annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = exon_tb.fetch(chr2, int(pos2) - 1, int(pos2))
    except Exception as inst:
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
   
    exon2 = [];
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            exon2.append(record[3])

    if len(exon2) == 0: exon2.append("---")
    exon2 = list(set(exon2))

    return [';'.join(gene1), ';'.join(gene2), ';'.join(exon1), ';'.join(exon2)]



def check_coding_info(chr, start, end, ref_coding_tb):

    # coding_score = {"coding": 4, "noncoding": 3, "3UTR": 2, "5UTR": 1, "intron": 0, "splicing": -1, "complex": -1}
    coding_score = {"coding": 11, "multi-coding": 10 , "coding_splicing": 9, "noncoding": 8, "noncoding_splicing": 7, 
                    "3UTR": 6, "3UTR_splicing": 5, "5UTR": 4, "5UTR_splicing": 3, "intron": 2, "complex": 1}

    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = ref_coding_tb.fetch(chr, int(start) - 1, int(end) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    coding_info = {}
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(start) + 1 <= int(record[2]) and int(end) - 1 >= int(record[1]) + 1: 
                coding_info[','.join([record[3], record[0], record[1], record[2], record[5]])] = record[4]


    gene2info = {}
    for key in sorted(coding_info):
        key_elm = key.split(',')
        gene = key_elm[0]
        gene2info[gene] = coding_info[key] if gene not in gene2info else gene2info[gene] + ';' + coding_info[key]

   
    info_bars = ""
    temp_type = "" 
    for gene in sorted(gene2info):
        info_bars = info_bars + ',' + gene + ';' + gene2info[gene]
        if len(gene2info[gene].split(';')) == 1:
            if temp_type == "": 
                temp_type = gene2info[gene]
            else:
                if coding_score[gene2info[gene]] > coding_score[temp_type]:
                    temp_type = gene2info[gene]
        elif len(gene2info[gene].split(';')) == 2 and gene2info[gene].find("intron") != -1:
            if gene2info[gene].find("coding") != -1: temp_type = "coding_splicing"
            elif gene2info[gene].find("noncoding") != -1: temp_type = "noncoding_splicing"
            elif gene2info[gene].find("3UTR") != -1: temp_type = "3UTR_splicing"
            elif gene2info[gene].find("5UTR") != -1: temp_type = "5UTR_splicing"
            else: temp_type = "complex"
        # elif len(list(set(gene2info[gene].split(';')))) == 2
        else:
            temp_type = "complex"
          
                  
    return temp_type + '\t' + info_bars.lstrip(',') 


def check_coding_info2(chr, start, end, ref_coding_tb):

    coding_score = {"C": 11, "C;I": 10, "I;C": 10, "I;C;I": 9, "C;I;C": 8, "I;C;I;C": 8, "C;I;C;I": 8, "I;C;I;C;I": 7, "C;I;C;I;C": 6, "I;C;I;C;I;C": 6, "C;I;C;I;C;I": 6,
                    "N": 5, "3": 4, "3;I": 3, "I;3": 3, "5": 2, "5;I": 1, "I;5": 1, "I": 0, "other": -1}

    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = ref_coding_tb.fetch(chr, int(start) - 1, int(end) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    coding_info = {}
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(start) + 1 <= int(record[2]) and int(end) - 1 >= int(record[1]) + 1:
                coding_info[','.join([record[3], record[0], record[1], record[2], record[5]])] = record[4][0].upper()


    gene2info = {}
    for key in sorted(coding_info):
        key_elm = key.split(',')
        gene = key_elm[0]
        gene2info[gene] = coding_info[key] if gene not in gene2info else gene2info[gene] + ';' + coding_info[key]


    info_bars = ""
    temp_type = "other"
    for gene in sorted(gene2info):
        info_bars = info_bars + ',' + gene + ';' + gene2info[gene]
        if gene2info[gene] in coding_score:
            if temp_type == "other":
                temp_type = gene2info[gene]
            else:
                if coding_score[gene2info[gene]] > coding_score[temp_type]: temp_type = gene2info[gene]

    return temp_type + '\t' + info_bars.lstrip(',')


def check_fusion_direction(chr1, pos1, dir1, chr2, pos2, dir2, ref_gene_tb, fusion_info_file):

    fusion_comb2info = {}
    with open(fusion_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            fusion_comb2info[F[0] + '\t' + F[1]] = F[2]

    potential_upstream_gene1 = []
    potential_downstream_gene1 = []
    potential_upstream_gene2 = []
    potential_downstream_gene2 = []

    ##########
    # check gene annotation for the side 1  
    tabixErrorFlag = 0
    try:
        records = ref_gene_tb.fetch(chr1, int(pos1) - 1, int(pos1) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(pos1) + 1 <= int(record[2]) and int(pos1) - 1 >= int(record[1]) + 1:
                if dir1 == record[5]: 
                    potential_upstream_gene1.append(record[3])
                else:
                    potential_downstream_gene1.append(record[3])

    ##########
    # check gene annotation for the side 2
    tabixErrorFlag = 0
    try:
        records = ref_gene_tb.fetch(chr2, int(pos2) - 1, int(pos2) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(pos2) + 1 <= int(record[2]) and int(pos2) - 1 >= int(record[1]) + 1:
                if dir2 == record[5]: 
                    potential_upstream_gene2.append(record[3])
                else:
                    potential_downstream_gene2.append(record[3])


    fusion_comb = []
    for u_gene in potential_upstream_gene1:
        for d_gene in potential_downstream_gene2:
            if u_gene != d_gene: fusion_comb.append(u_gene + ';' + d_gene)
 
    for u_gene in potential_upstream_gene2:
        for d_gene in potential_downstream_gene1:
            if u_gene != d_gene: fusion_comb.append(u_gene + ';' + d_gene)


    fusion_comb = list(set(fusion_comb))

    fusion_infos = [] 
    for comb in fusion_comb:
        g1, g2 = sorted(comb.split(';'))
        if g1 + '\t' + g2 in fusion_comb2info: 
            for db in fusion_comb2info[g1 + '\t' + g2].split(';'):
                fusion_infos.append(db)

    fusion_info = ';'.join(sorted(list(set(fusion_infos)))) if len(fusion_infos) > 0 else "---"

    return ','.join(fusion_comb) + '\t' + fusion_info if len(fusion_comb) > 0 else "---\t---"



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


