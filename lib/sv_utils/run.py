#! /usr/bin/env python

import os, subprocess
import pysam
import utils

def count_main(args):

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
                                                args.max_size_thres, ref_junc_tb, ens_junc_tb, grch2ucsc)

            if len(sv_good_list) > 0:
                for i in range(0, len(sv_good_list)):
                    type2count[sv_good_list[i][7]] = type2count[sv_good_list[i][7]] + 1

            total = type2count["deletion"] + type2count["tandem_duplication"] + type2count["inversion"] + type2count["translocation"]
            print >> hout, sample + '\t' + tumor_type + '\t' + str(type2count["deletion"]) + '\t' + str(type2count["tandem_duplication"]) + '\t' + \
                                        str(type2count["inversion"]) + '\t' + str(type2count["translocation"]) + '\t' + str(total) 
 
    hout.close()


def gene_summary_main(args):

    hout = open(args.output, 'w')

    # read cancer gene
    gene2info = {}
    with open(args.cancer_gene_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene2info[F[0]] = '\t'.join(F[1:])


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

    tumor_type_list = {} 
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres,
                                                args.normal_freq_thres, args.normal_depth_thres, args.inversion_size_thres,
                                                args.max_size_thres, ref_junc_tb, ens_junc_tb, grch2ucsc)

            if len(sv_good_list) > 0:
                for i in range(0, len(sv_good_list)):
                    var_type = sv_good_list[i][7]
                    genes1 = sv_good_list[i][8].split(';') if sv_good_list[i][8] != "---" else []
                    genes2 = sv_good_list[i][9].split(';') if sv_good_list[i][9] != "---" else []
                    for gene in list(set(genes1 + genes2)):
                        if gene in gene2type_sample:
                            gene2type_sample[gene] = gene2type_sample[gene] + [(tumor_type, sample, var_type)]
                        else:
                            gene2type_sample[gene] = [(tumor_type, sample, var_type)]


    print >> hout, '\t'.join(["gene", "all_type_total", "all_type_each"] + \
                reduce(lambda x, y: x + y, [[x + "_total", x + "_each"] for x in sorted(tumor_type_list.keys())]) + ["Lawrence et al", "CGC"])

    for gene in sorted(gene2type_sample):
        tumor_sample_var = list(set(gene2type_sample[gene]))
        tumor2count= {}
        total_count = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0}
        for tumor_type in tumor_type_list.keys(): tumor2count[tumor_type] = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0} 
        for tumor_type, sample, var in tumor_sample_var:
            tumor2count[tumor_type][var] = tumor2count[tumor_type][var] + 1
            total_count[var] = total_count[var] + 1

        count_bar = ""
        for tumor_type in sorted(tumor2count):
            count_bar = count_bar + '\t' + str(tumor2count[tumor_type]["deletion"] + tumor2count[tumor_type]["tandem_duplication"] + tumor2count[tumor_type]["inversion"] + tumor2count[tumor_type]["translocation"]) + '\t' + \
                        str(tumor2count[tumor_type]["deletion"]) + ',' + str(tumor2count[tumor_type]["tandem_duplication"]) + ',' + str(tumor2count[tumor_type]["inversion"]) + ',' + str(tumor2count[tumor_type]["translocation"])

        info = gene2info[gene] if gene in gene2info else "---" + '\t' + "---"
        print >> hout, gene + '\t' + str(total_count["deletion"] + total_count["tandem_duplication"] + total_count["inversion"] + total_count["translocation"]) + '\t' + \
                str(total_count["deletion"]) + ',' + str(total_count["tandem_duplication"]) + ',' + str(total_count["inversion"]) + ',' + str(total_count["translocation"]) +  count_bar + '\t' + info

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
                                  args.max_size_thres, ref_junc_tb, ens_junc_tb, grch2ucsc)

    if len(sv_good_list) > 0:
        for i in range(0, len(sv_good_list)):
            print >> hout, '\t'.join(sv_good_list[i])

    hout.close()


def concentrate_main(args):

    hout = open(args.output + ".tmp.sv_list", 'w')

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

    tumor_type_list = {}
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres,
                                                args.normal_freq_thres, args.normal_depth_thres, args.inversion_size_thres,
                                                args.max_size_thres, ref_junc_tb, ens_junc_tb, grch2ucsc)

            if len(sv_good_list) > 0:
                for i in range(0, len(sv_good_list)):
                    if sv_good_list[i][7] in ["deletion", "tandem_duplication"]:
                        print >> hout, sample + '\t' + tumor_type + '\t' + '\t'.join(sv_good_list[i])

    hout.close()

    hout = open(args.output + ".tmp.sv_list.sort", 'w')
    subprocess.call(["sort", "-k3,3", "-k4,4n", "-k6,6", "-k7,7n", args.output + ".tmp.sv_list"], stdout = hout)
    hout.close()


    margin = args.set_margin
    count_thres = args.set_count 
    temp_keys = []
    temp_ind = -1
    passed_keys = []

    hout = open(args.output, 'w')

    with open(args.output + ".tmp.sv_list.sort", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[3] == "178098812":
                pass

            temp_keys.append(F)
            if len(temp_keys) < count_thres: continue

            if temp_keys[0][2] == temp_keys[count_thres - 1][2] and int(temp_keys[count_thres - 1][3]) - int(temp_keys[0][3]) <= margin:
                # current frame satisties the concentrated condition
                if temp_ind == -1: 
                    passed_keys = passed_keys + temp_keys
                else:
                    passed_keys = passed_keys + temp_keys[temp_ind:]
                temp_ind = count_thres - 1
            else:
                # current frame does not satisfy the concentrated condition
                if temp_ind == -1:
                    # flush
                    for i in range(0, len(passed_keys)):
                        print >> hout, '\t'.join(passed_keys[i])
                    passed_keys = []
                    temp_ind = -1
                else:
                    temp_ind = temp_ind - 1

            temp_keys.pop(0)

        if len(passed_keys) > 0:
            for i in range(0, len(passed_keys)):
                print '\t'.join(passed_keys[i])

    hout.close()

    subprocess.call(["rm", "-rf", args.output + ".tmp.sv_list"])
    subprocess.call(["rm", "-rf", args.output + ".tmp.sv_list.sort"])

