#! /usr/bin/env python

import sys, pkg_resources


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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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


def check_fusion_direction(chr1, pos1, dir1, chr2, pos2, dir2, ref_gene_tb, known_fusion_info_file):

    # known_fusion_info_file = pkg_resources.resource_filename("sv_utils", "data/fusion_db.txt")

    known_fusion_comb2info = {}
    with open(known_fusion_info_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            known_fusion_comb2info[F[0] + '\t' + F[1]] = F[2]

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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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
        # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
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

    known_fusion_comb = []
    known_fusion_source = [] 
    for comb in fusion_comb:
        g1, g2 = sorted(comb.split(';'))
        if g1 + '\t' + g2 in known_fusion_comb2info:
            known_fusion_comb.append(g1 + ';' + g2)
            known_fusion_source.append(known_fusion_comb2info[g1 + '\t' + g2])

    if len(known_fusion_comb) > 0:
        return ','.join(known_fusion_comb) + '\t' + ','.join(known_fusion_source)
    else:
        return "---\t---"



