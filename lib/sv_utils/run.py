#! /usr/bin/env python

import os
import pysam
import utils

def summary_main(args):

    hout = open(args.output, 'w')
    print >> hout, '\t'.join(["sample", "type", "deletion", "tandem_duplication", "inversion", "translocation", "total"])

    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]


    with open(args.result_list, 'r') as hin:

        for line in hin:
            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            type2count = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0}

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres,
                                                args.normal_freq_thres, args.normal_depth_thres, args.inversion_size_thres,
                                                ref_junc_tb, ens_junc_tb, grch2ucsc)

            if len(sv_good_list) > 0:
                for i in range(0, len(sv_good_list)):
                    type2count[sv_good_list[i][7]] = type2count[sv_good_list[i][7]] + 1

            total = type2count["deletion"] + type2count["tandem_duplication"] + type2count["inversion"] + type2count["translocation"]
            print >> hout, sample + '\t' + tumor_type + '\t' + str(type2count["deletion"]) + '\t' + str(type2count["tandem_duplication"]) + '\t' + \
                                        str(type2count["inversion"]) + '\t' + str(type2count["translocation"]) + '\t' + str(total) 
 
    hout.close()


def filter_main(args):

    hout = open(args.output, 'w')

    if not os.path.exists(args.result_file):
        raise ValueError("file not exists: " + args.result_file)

    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    sv_good_list = utils.filter_sv_list(args.result_file, args.fisher_thres, args.tumor_freq_thres,
                                  args.normal_freq_thres, args.normal_depth_thres, args.inversion_size_thres,
                                  ref_junc_tb, ens_junc_tb, grch2ucsc)

    if len(sv_good_list) > 0:
        for i in range(0, len(sv_good_list)):
            print >> hout, '\t'.join(sv_good_list[i])

    hout.close()


