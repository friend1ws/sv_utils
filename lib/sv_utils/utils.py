#! /usr/bin/env python

import sys, subprocess
import my_seq

def filter_sv_list(result_file, fisher_thres, tumor_freq_thres, normal_freq_thres, normal_depth_thres, 
                    inversion_size_thres, max_size_thres, within_exon, ref_exon_tb, ens_exon_tb, ref_junc_tb, ens_junc_tb, 
                    simple_repeat_tb, grch2ucsc, control_tb, control_num_thres, normal_mode = False):

    good_list = []
    with open(result_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[0] == "MT" or F[3] == "MT": continue
            if F[0] == "hs37d5" or F[3] == "hs37d5": continue
            if F[0].startswith("GL0") or F[3].startswith("GL0"): continue

            if F[7] == "inversion" and abs(int(F[1]) - int(F[4])) < int(inversion_size_thres): continue
            if max_size_thres is not None:
                if F[7] == "translocation": continue
                if abs(int(F[1]) - int(F[4])) > int(max_size_thres): continue
            if float(F[14]) < float(tumor_freq_thres): continue

            if normal_mode == False:
                if int(F[15]) + int(F[16]) < int(normal_depth_thres): continue
                if float(F[17]) > float(normal_freq_thres): continue
                if float(F[18]) < float(fisher_thres): continue

            if within_exon:
                if F[7] == "translocation": continue
                chr_ucsc = grch2ucsc[F[0]] if F[0] in grch2ucsc else F[0]
                if not in_exon_check(chr_ucsc, F[1], F[4], ref_exon_tb, ens_exon_tb): continue
    
            if control_tb is not None:
                if control_check(F[0], F[1], F[2], F[3], F[4], F[5], F[6], control_tb, control_num_thres): 
                    continue

            if F[7] in ["deletion", "tandem_duplication"] and simple_repeat_tb is not None:
                chr_ucsc = grch2ucsc[F[0]] if F[0] in grch2ucsc else F[0]
                if simple_repeat_check(chr_ucsc, F[1], F[4], simple_repeat_tb):
                    continue

            if F[7] == "deletion":
                chr_ucsc = grch2ucsc[F[0]] if F[0] in grch2ucsc else F[0]
                if junction_check(chr_ucsc, F[1], F[4], ref_junc_tb, ens_junc_tb): continue

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
def junction_check(chr, start, end, ref_junc_tb, ens_junc_tb, margin = 1):

    junc_flag = False
 
    # check junction annotation for refGene 
    tabixErrorFlag = 0
    try:
        records = ref_junc_tb.fetch(chr, int(start) - 30, int(start) + 30)
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
        records = ens_junc_tb.fetch(chr, int(start) - 30, int(start) + 30)
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
def distance_to_closest(chr1, pos1, chr2, pos2, ref_exon_tb, search_max = 500000):

    cur_dist = search_max
    target_gene = []

    # check junction annotation for refGene for the first break point
    tabixErrorFlag = 0
    try:
        search_start = max(0, int(pos1) - search_max)
        records = ref_exon_tb.fetch(chr1, search_start, int(pos1) + search_max)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1

    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
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
        records = ref_exon_tb.fetch(chr2, search_start, int(pos2) + search_max)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag = 1
    
    if tabixErrorFlag == 0:
        for record_line in records:
            record = record_line.split('\t')
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
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(0, len(header)):
            if header[i] == "Ref": ref_ind = i
            if header[i] == "Alt": alt_ind = i
            if header[i] == "readPairNum_tumor": tum_ref_ind = i
            if header[i] == "variantPairNum_tumor": tum_var_ind = i
            if header[i] == "readPairNum_normal": nor_ref_ind = i
            if header[i] == "variantPairNum_normal": nor_var_ind = i
            if header[i] == "P-value(fhsher_realignment)": fisher_ind = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if len(F[ref_ind]) >= 10 or len(F[alt_ind]) >= 10:

                bed_key = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[2]
                read_info = F[tum_ref_ind] + '\t' + F[tum_var_ind] + '\t' + \
                                str(round(float(F[tum_var_ind]) / (float(F[tum_ref_ind]) + float(F[tum_var_ind])), 3)) + '\t' + \
                                F[nor_ref_ind] + '\t' + F[nor_var_ind] + '\t' + \
                                str(round(float(F[nor_var_ind]) / (float(F[nor_ref_ind]) + float(F[nor_var_ind])), 3)) + '\t' + F[fisher_ind]

                var_info = ""
                # deletion
                if len(F[ref_ind]) >= 10:
                    var_info = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + "+" + '\t' + \
                                   F[0] + '\t' + str(int(F[1]) + len(F[ref_ind])) + '\t' + "-" + '\t' + "---" + '\t' + "deletion"
                    # gene_annotation = get_gene_annotation(F[0],  str(int(F[1]) - 1), F[0], str(int(F[1]) + len(F[ref_ind])), gene_tb, exon_tb)
                    gene_annotation = "---" + '\t' + "---" + '\t' + "---" + '\t' + "---"
                    print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info
 
                # tandem_duplication
                elif len(F[alt_ind]) >= 10:
                    # tandem_duplication check
                    flanking_seq_1 = my_seq.get_seq(reference, F[0], int(F[1]) + 1, int(F[1]) + len(F[alt_ind]))
                    flanking_seq_1_match = my_seq.exact_alignment(F[alt_ind], flanking_seq_1)
                    flanking_seq_2 = my_seq.get_seq(reference, F[0], int(F[1]) - len(F[alt_ind]) + 1, int(F[1])) 
                    flanking_seq_2_match = my_seq.exact_alignment(F[alt_ind], flanking_seq_2)
                    # print '\t'.join(F[0:4])
                    # print F[alt_ind] + '\t' + flanking_seq_1 + '\t' + str(flanking_seq_1_match)
                    # print F[alt_ind] + '\t' + flanking_seq_2 + '\t' + str(flanking_seq_2_match)
    
                    if flanking_seq_1_match == len(F[alt_ind]) or flanking_seq_2_match == len(F[alt_ind]):
                        var_info = F[0] + '\t' +  str(int(F[1]) + 1) + '\t' + "-" + '\t' + \
                                       F[0] + '\t' + str(int(F[1]) + len(F[alt_ind])) + '\t' + "+" + '\t' + "---" + '\t' + "tandem_duplication"
                        # gene_annotation = get_gene_annotation(F[0],  str(int(F[1]) + 1), F[0], str(int(F[1]) + len(F[alt_ind])), gene_tb, exon_tb)
                        gene_annotation = "---" + '\t' + "---" + '\t' + "---" + '\t' + "---"
                        print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info
 
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

    coding_score = {"coding": 4, "noncoding": 3, "3UTR": 2, "5UTR": 1, "intron": 0, "spliing": -1, "complex": -1}

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
        elif len(gene2info[gene].split(';')) == 2 and gene2info[gene].find("intron"):
            temp_type = "splicing"
        else:
            temp_type = "complex"
          
                  
    return temp_type + '\t' + info_bars.lstrip(',') 



