#! /usr/bin/env python

import sys, subprocess
from utils import header_check
from simple_repeat import make_simple_repeat_info
from header_info import *
from annot_utils.junction import *
# from annot_utils.simple_repeat import *
import pysam

# def filter_sv_list(args, ref_exon_tb, ens_exon_tb, ref_junc_tb, ens_junc_tb, 
#                     simple_repeat_tb, grch2ucsc, control_tb):

def filter_sv_list(args):


    if args.remove_rna_junction:
        refseq_junc_info = args.output_file + ".refseq.junc.bed.gz"
        gencode_junc_info = args.output_file + ".gencode.junc.bed.gz"

        make_junc_info(refseq_junc_info, "refseq", args.genome_id, args.grc, False)
        make_junc_info(gencode_junc_info, "gencode", args.genome_id, args.grc, False)
        refseq_junc_tb = pysam.TabixFile(refseq_junc_info)
        gencode_junc_tb = pysam.TabixFile(gencode_junc_info)

    if args.simple_repeat_file is not None:
        simple_repeat_info = args.output_file + ".simple_repeat.bed.gz"
        make_simple_repeat_info(args.simple_repeat_file, simple_repeat_info, args.genome_id, args.grc)
        simple_repeat_tb = pysam.TabixFile(simple_repeat_info)

    control_tb = pysam.TabixFile(args.pooled_control_file) if args.pooled_control_file is not None else None


    good_list = []
    with open(args.input_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"):
                good_list.append(line.rstrip('\n'))
                continue

            if header_check(line.rstrip('\n')):
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


            """
            if args.within_exon:
                if F[header_info.variant_type] == "translocation": continue
                chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if not in_exon_check(chr_ucsc, F[header_info.pos_1], F[header_info.pos_2], ref_exon_tb, ens_exon_tb): continue
            """
    
            if control_tb is not None:
                if control_check(F[header_info.chr_1], F[header_info.pos_1], F[header_info.dir_1], \
                                 F[header_info.chr_2], F[header_info.pos_2], F[header_info.dir_2], \
                                 F[header_info.inserted_seq], control_tb, args.pooled_control_num_thres): 
                    continue

            if F[header_info.variant_type] in ["deletion", "tandem_duplication"] and args.simple_repeat_file is not None:
                # chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if simple_repeat_check(F[header_info.chr_1], F[header_info.pos_1], F[header_info.pos_2], simple_repeat_tb): continue

            if F[header_info.variant_type] == "deletion" and args.remove_rna_junction:
                # chr_ucsc = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1] 
                if junction_check(F[header_info.chr_1], F[header_info.pos_1], F[header_info.pos_2], refseq_junc_tb, gencode_junc_tb): continue

            good_list.append(F)

    
    if args.remove_rna_junction:
        subprocess.check_call(["rm", "-rf", refseq_junc_info])
        subprocess.check_call(["rm", "-rf", refseq_junc_info + ".tbi"])
        subprocess.check_call(["rm", "-rf", gencode_junc_info])
        subprocess.check_call(["rm", "-rf", gencode_junc_info + ".tbi"])

    if args.simple_repeat_file is not None:
        subprocess.check_call(["rm", "-rf", simple_repeat_info])
        subprocess.check_call(["rm", "-rf", simple_repeat_info + ".tbi"])

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


