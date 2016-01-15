#! /usr/bin/env python

def filter_sv_list(result_file, fisher_thres, tumor_freq_thres, normal_freq_thres, normal_depth_thres, 
                    inversion_size_thres, max_size_thres, ref_junc_tb, ens_junc_tb, grch2ucsc):

    good_list = []
    with open(result_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[7] == "inversion" and abs(int(F[1]) - int(F[4])) < int(inversion_size_thres): continue
            if max_size_thres is not None:
                if F[7] == "translocation": continue
                if abs(int(F[1]) - int(F[4])) > int(max_size_thres): continue
            if float(F[14]) < float(tumor_freq_thres): continue
            if int(F[15]) + int(F[16]) < int(normal_depth_thres): continue
            if float(F[17]) > float(normal_freq_thres): continue
            if float(F[18]) < float(fisher_thres): continue

            if F[7] == "deletion":
                chr_ucsc = grch2ucsc[F[0]] if F[0] in grch2ucsc else F[0]
                if junction_check(chr_ucsc, F[1], F[4], ref_junc_tb, ens_junc_tb): continue

            good_list.append(F)

    return good_list


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

    # check junction annotation for refGene 
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
 


