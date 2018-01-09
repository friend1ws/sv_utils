#! /usr/bin/env python


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


