#! /usr/bin/env python

import sys

def filter_sv_list(result_file, fisher_thres, tumor_freq_thres, normal_freq_thres, normal_depth_thres, 
                    inversion_size_thres, max_size_thres, within_exon, ref_exon_tb, ens_exon_tb, ref_junc_tb, ens_junc_tb, grch2ucsc, 
                    control_tb, control_num_thres, normal_mode = False):

    good_list = []
    with open(result_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

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
        records = ref_exon_tb.fetch(chr1, int(pos1) - search_max, int(pos1) + search_max)
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
        records = ref_exon_tb.fetch(chr2, int(pos2) - search_max, int(pos2) + search_max)
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



