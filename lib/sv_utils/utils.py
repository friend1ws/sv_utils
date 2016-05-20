#! /usr/bin/env python

import sys, subprocess
import my_seq
from header_info import *

def header_check(line):

    if "Chr_1" + '\t' + "Pos_1" + '\t' + "Dir_1" + '\t' + "Chr_2" + '\t' + "Pos_2" + '\t' + "Dir_2" + '\t' + "Inserted_Seq" + '\t' + "Variant_Type" in line:
        return True
    else:
        return False

def filter_sv_list(args, ref_exon_tb, ens_exon_tb, ref_junc_tb, ens_junc_tb, 
                    simple_repeat_tb, grch2ucsc, control_tb):

    good_list = []
    with open(args.result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"):
                good_list.append(line.rstrip('\n'))
                continue

            if utils.header_check(line.rstrip('\n')):
                header_info.read(line.rstrip('\n'))
                F = line.rstrip('\n').split('\t')
                good_list.append(F)
                continue

            F = line.rstrip('\n').split('\t')

            if F[header_info.chr_1] == "MT" or F[header_info.chr_2] == "MT": continue
            if [header_info.chr_1] == "hs37d5" or F[header_info.chr_2] == "hs37d5": continue
            if F[header_info.chr_1].startswith("GL0") or F[header_info.chr_2].startswith("GL0"): continue

            if F[header_info.variant_type] == "inversion" and abs(int(F[header_info.pos_1]) - int(F[header_info.pos_2])) < int(args.inversion_size_thres): continue
            if args.max_variant_size is not None:
                if F[header_info.variant_type] == "translocation": continue
                if abs(int(F[header_info.pos_1]) - int(F[header_info.pos_2])) > int(args.args.max_variant_size): continue
            if float(F[header_info.tumor_vaf]) < float(args.min_tumor_allele_freq): continue

            if F[header_info.num_control_var_read_pair] != "---":
                if int(F[header_info.num_control_var_read_pair]) > int(args.max_control_variant_read_pair): continue

            if F[header_info.num_control_ref_read_pair] != "---" and F[header_info.num_control_var_read_pair] != "---":    
                if int(F[header_info.num_control_ref_read_pair]) + int(F[header_info.num_control_var_read_pair]) < int(args.control_depth_thres): continue

            if F[header_info.control_vaf] != "---":
                if float(F[header_info.control_vaf]) > float(args.max_control_allele_freq): continue

            if F[header_info.minus_log_fisher_p_value] != "---":
                if float(F[header_info.minus_log_fisher_p_value]) < float(args.max_minus_log_fisher_pvalue): continue

            if int(F[header_info.max_over_hang_1]) < int(args.min_overhang_size) or int(F[header_info.max_over_hang_2]) < int(args.min_overhang_size): continue


            if args.within_exon:
                if F[header_info.variant_type] == "translocation": continue
                chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if not in_exon_check(chr_ucsc, F[header_info.pos_1], F[header_info.pos_2], ref_exon_tb, ens_exon_tb): continue
    
            if control_tb is not None:
                if control_check(F[header_info.chr_1], F[header_info.pos_1], F[header_info.dir_1], \
                                 F[header_info.chr_2], F[header_info.pos_2], F[header_info.dir_2], \
                                 F[header_info.inserted_seq], control_tb, args.pooled_control_num_thres): 
                    continue

            if F[header_info.variant_type] in ["deletion", "tandem_duplication"] and simple_repeat_tb is not None:
                chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if simple_repeat_check(chr_ucsc, F[header_info.pos_1], F[header_info.pos_2], simple_repeat_tb):
                    continue

            if F[header_info.variant_type] == "deletion":
                chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if junction_check(chr_ucsc, F[header_info.pos_1], F[header_info.pos_2], ref_junc_tb, ens_junc_tb): continue

            good_list.append(F)

    return good_list


def control_check(chr1, pos1, dir1, chr2, pos2, dir2, inseq, control_tb, control_num_thres):

    inseq_size = (0 if inseq == "---" else len(inseq))
    control_panel_check_margin = 50 # should this be changable?
    control_flag = 0

    ####################
    # get the records for control junction data for the current position
    tabixErrorFlag = 0
    try:
        records = control_tb.fetch(chr1, int(pos1) - control_panel_check_margin, int(pos1) + control_panel_check_margin)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorMsg = str(inst.args)
        tabixErrorFlag = 1
    
    ####################

       ####################
    # for each record in control junction extracted, check the consistency with the current junction
    if tabixErrorFlag == 0:

        match_record_count = 0
        for record_line in records:
            record = record_line.split('\t')

            if chr1 == record[0] and chr2 == record[3] and dir1 == record[8] and dir2 == record[9]:
    
                record_inseq_size = (0 if record[7] == "---" else len(record[7]))

                # detailed check on the junction position considering inserted sequences
                if dir1 == "+":
                    expectedDiffSize = (int(pos1) - int(record[2])) + (inseq_size - record_inseq_size)
                    if (dir2 == "+" and int(pos2) == int(record[5]) - int(expectedDiffSize)) or (dir2 == "-" and int(pos2) == int(record[5]) + int(expectedDiffSize)):
                        match_record_count = match_record_count + 1
                else:
                    expectedDiffSize = (int(pos1) - int(record[2])) + (record_inseq_size - inseq_size)
                    if (dir2 == "+" and int(pos2) == int(record[5]) + int(expectedDiffSize)) or (dir2 == "-" and int(pos2) == int(record[5]) - int(expectedDiffSize)):
                        match_record_count = match_record_count + 1

        return match_record_count >= control_num_thres

    else:
        
        return False



def in_exon_check(chr, start, end, ref_exon_tb, ens_exon_tb, margin = 5):

    in_exon_flag = False

    # check exon annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_exon_tb.fetch(chr, int(start) - margin * 2, int(start) + margin * 2)
    except Exception as inst:
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            exon_start = int(record[1])
            exon_end = int(record[2])
            if exon_start + 1 - margin <= int(start) and int(end) <= exon_end + margin: in_exon_flag = True


    # check exon annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ens_exon_tb.fetch(chr, int(start) - margin * 2, int(start) + margin * 2)
    except Exception as inst:
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            exon_start = int(record[1]) 
            exon_end = int(record[2]) 
            if exon_start + 1 - margin <= int(start) and int(end) <= exon_end + margin: in_exon_flag = True

    return in_exon_flag


# check whether potential rna contamination or not
def junction_check(chr, start, end, ref_junc_tb, ens_junc_tb, margin = 2):

    junc_flag = False
 
    # check junction annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_junc_tb.fetch(chr, int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: junc_flag = True

    # check junction annotation for ensGene 
    tabixErrorFlag = 0
    try:
        records = ens_junc_tb.fetch(chr, int(start) - 10, int(start) + 10)
    except Exception as inst:
        tabixErrorFlag = 1
        
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            start_diff = int(record[1]) - int(start)
            end_diff = int(record[2]) + 1 - int(end)
            if abs(start_diff - end_diff) <= margin: junc_flag = True

    return junc_flag
 



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



# def make_mut_db(input_file, output_file_prefix, reference, gene_tb, exon_tb):
def make_mut_db(input_file, output_file_prefix, reference):

    hout = open(output_file_prefix + ".bed", 'w')
    with open(input_file, 'r') as hin:

        ref_ind = -1
        alt_ind = -1
        tum_ref_ind = -1
        nor_ref_ind = -1
        tum_var_ind = -1 
        nor_var_ind = -1
        fisher_ind = -1

        header_line = "#"
        while header_line.startswith("#"):
            header_line = hin.readline().rstrip('\n')

        header = header_line.split('\t')

        for i in range(0, len(header)):
            if header[i] == "Ref": ref_ind = i
            if header[i] == "Alt": alt_ind = i
            if header[i] == "readPairNum_tumor": tum_ref_ind = i
            if header[i] == "variantPairNum_tumor": tum_var_ind = i
            if header[i] == "readPairNum_normal": nor_ref_ind = i
            if header[i] == "variantPairNum_normal": nor_var_ind = i
            if header[i] == "P-value(fisher_realignment)": fisher_ind = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if len(F[ref_ind]) >= 10 or len(F[alt_ind]) >= 10:

                if F[nor_ref_ind] == "---": continue
                if F[nor_var_ind] == "---": continue
                if F[tum_ref_ind] == "---": continue
                if F[tum_var_ind] == "---": continue

                bed_key = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[2]
                read_info = F[tum_ref_ind] + '\t' + F[tum_var_ind] + '\t' + \
                                str(round(float(F[tum_var_ind]) / (float(F[tum_ref_ind]) + float(F[tum_var_ind])), 3)) + '\t' + \
                                F[nor_ref_ind] + '\t' + F[nor_var_ind] + '\t' + \
                                str(round(float(F[nor_var_ind]) / (float(F[nor_ref_ind]) + float(F[nor_var_ind])), 3)) + '\t' + F[fisher_ind]

                var_info = ""
                gene_annotation = "---" + '\t' + "---" + '\t' + "---" + '\t' + "---"
                other_info = "---" + '\t' + "0" + '\t' + "---" + '\t' + "---"
                # deletion
                if len(F[ref_ind]) >= 10:
                    var_info = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + "+" + '\t' + \
                                   F[0] + '\t' + str(int(F[1]) + len(F[ref_ind])) + '\t' + "-" + '\t' + "---" + '\t' + "deletion"
                    print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info + '\t' + other_info
 
                # tandem_duplication
                elif len(F[alt_ind]) >= 10:
                    # tandem_duplication check
                    flanking_seq_1 = my_seq.get_seq(reference, F[0], int(F[1]) + 1, int(F[1]) + len(F[alt_ind]))
                    flanking_seq_1_match, coord1 = my_seq.exact_alignment(F[alt_ind], flanking_seq_1)
                    flanking_seq_2 = my_seq.get_seq(reference, F[0], int(F[1]) - len(F[alt_ind]) + 1, int(F[1])) 
                    flanking_seq_2_match, coord2 = my_seq.exact_alignment(F[alt_ind], flanking_seq_2)
    
                    if flanking_seq_1_match == len(F[alt_ind]) or flanking_seq_2_match == len(F[alt_ind]):
                        var_info = F[0] + '\t' +  str(int(F[1]) + 1) + '\t' + "-" + '\t' + \
                                       F[0] + '\t' + str(int(F[1]) + len(F[alt_ind])) + '\t' + "+" + '\t' + "---" + '\t' + "tandem_duplication"
                        print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info + '\t' + other_info
 
    hout.close()

    subprocess.call(["bgzip", "-f", output_file_prefix + ".bed"])
    subprocess.call(["tabix", "-p", "bed", output_file_prefix + ".bed.gz"])



def simple_repeat_check(chr, start, end, simple_repeat_tb):

    # check junction annotation for refGene for the first break point
    tabixErrorFlag = 0
    try:
        records = simple_repeat_tb.fetch(chr, int(start) - 1, int(end) + 1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
            if int(record[1]) <= int(start) + 1 and int(end) - 1 <= int(record[2]):
                return True

    return False



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


