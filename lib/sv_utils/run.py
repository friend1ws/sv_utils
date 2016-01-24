#! /usr/bin/env python

import sys, os, subprocess
import pysam
import utils

def count_main(args):


    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    
    hout = open(args.output, 'w')
    if args.inseq == True:
        print >> hout, '\t'.join(["sample", "type", "deletion_nonseq", "deletion_inseq", "tandem_duplication_nonseq", "tandem_duplication_inseq",
                                    "inversion_nonseq", "inversion_inseq", "translocation_nonseq", "translocation_inseq", "total_nonseq", "total_inseq"])
    else:
        print >> hout, '\t'.join(["sample", "type", "deletion", "tandem_duplication", "inversion", "translocation", "total"])


    with open(args.result_list, 'r') as hin:

        for line in hin:
            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres, args.normal_freq_thres,
                                                 args.normal_depth_thres, args.inversion_size_thres, args.max_size_thres,
                                                 args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_bed, ens_junc_bed, grch2ucsc, 
                                                 control_tb, args.control_num_thres, False)
            
            if args.inseq == True:
                type2count = {"deletion_nonseq": 0, "deletion_inseq": 0, "tandem_duplication_nonseq": 0, "tandem_duplication_inseq": 0,
                                "inversion_nonseq": 0, "inversion_inseq": 0, "translocation_nonseq": 0, "translocation_inseq": 0}

                if len(sv_good_list) > 0:
                    for i in range(0, len(sv_good_list)):
                        if sv_good_list[i][6] != "---":
                            type2count[sv_good_list[i][7] + "_inseq"] = type2count[sv_good_list[i][7] + "_inseq"] + 1
                        else:
                            type2count[sv_good_list[i][7] + "_nonseq"] = type2count[sv_good_list[i][7] + "_nonseq"] + 1

                total_nonseq = type2count["deletion_nonseq"] + type2count["tandem_duplication_nonseq"] + \
                                    type2count["inversion_nonseq"] + type2count["translocation_nonseq"]
                total_inseq = type2count["deletion_inseq"] + type2count["tandem_duplication_inseq"] + \
                                    type2count["inversion_inseq"] + type2count["translocation_inseq"]

                print >> hout, sample + '\t' + tumor_type + '\t' + str(type2count["deletion_nonseq"]) + '\t' + str(type2count["deletion_inseq"]) + '\t' + \
                                    str(type2count["tandem_duplication_nonseq"]) + '\t' + str(type2count["tandem_duplication_inseq"]) + '\t' + \
                                    str(type2count["inversion_nonseq"]) + '\t' + str(type2count["inversion_inseq"]) + '\t' + \
                                    str(type2count["translocation_nonseq"]) + '\t' + str(type2count["translocation_inseq"]) + '\t' + \
                                    str(total_nonseq) + '\t' + str(total_inseq)

            else:
                type2count = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0}

                if len(sv_good_list) > 0:
                    for i in range(0, len(sv_good_list)):
                        type2count[sv_good_list[i][7]] = type2count[sv_good_list[i][7]] + 1

                total = type2count["deletion"] + type2count["tandem_duplication"] + type2count["inversion"] + type2count["translocation"]
                print >> hout, sample + '\t' + tumor_type + '\t' + str(type2count["deletion"]) + '\t' + str(type2count["tandem_duplication"]) + '\t' + \
                                            str(type2count["inversion"]) + '\t' + str(type2count["translocation"]) + '\t' + str(total) 
 
    hout.close()


def gene_summary_main(args):


    # read cancer gene
    gene2info = {}
    with open(args.cancer_gene_list, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            gene2info[F[0]] = '\t'.join(F[1:])


    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    
    hout = open(args.output, 'w')

    tumor_type_list = {} 
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            print >> sys.stderr, "reading: " + sample + ' ' + tumor_type

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres, args.normal_freq_thres,
                                                args.normal_depth_thres, args.inversion_size_thres, args.max_size_thres,
                                                args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_bed, ens_junc_bed, grch2ucsc, 
                                                control_tb, args.control_num_thres, False)

            if len(sv_good_list) > 0:
                for i in range(0, len(sv_good_list)):
                    var_type = sv_good_list[i][7]
                    genes1 = sv_good_list[i][8].split(';') if sv_good_list[i][8] != "---" else []
                    genes2 = sv_good_list[i][9].split(';') if sv_good_list[i][9] != "---" else []
    
                    # check in-frame or not
                    is_inframe = False
                    if var_type == "deletion":
                        var_size = int(sv_good_list[i][4]) - int(sv_good_list[i][1]) - 1
                        var_size = var_size - (0 if sv_good_list[i][6] == "---" else len(sv_good_list[i][6]))
                        if var_size % 3 == 0: is_inframe = True
                    elif var_type == "tandem_duplication":
                        var_size = int(sv_good_list[i][4]) - int(sv_good_list[i][1]) + 1
                        var_size = var_size + (0 if sv_good_list[i][6] == "---" else len(sv_good_list[i][6]))
                        if var_size % 3 == 0: is_inframe = True

                    for gene in list(set(genes1 + genes2)):
                        if gene in gene2type_sample:
                            gene2type_sample[gene] = gene2type_sample[gene] + [(tumor_type, sample, var_type, is_inframe)]
                        else:
                            gene2type_sample[gene] = [(tumor_type, sample, var_type, is_inframe)]


    if args.inframe_info == True:
        print >> hout, '\t'.join(["gene", "all_type_total", "all_type_each", "all_type_total_inframe", "all_type_each_inframe"] + \
                    reduce(lambda x, y: x + y, [[x + "_total", x + "_each", x + "_total_inframe", x + "_each_inframe"] for x in sorted(tumor_type_list.keys())]) + \
                    ["Lawrence et al", "CGC"])
    else:
        print >> hout, '\t'.join(["gene", "all_type_total", "all_type_each"] + \
                    reduce(lambda x, y: x + y, [[x + "_total", x + "_each"] for x in sorted(tumor_type_list.keys())]) + ["Lawrence et al", "CGC", "kinbase"])

    for gene in sorted(gene2type_sample):
        tumor_sample_var = list(set(gene2type_sample[gene]))
        tumor2count= {}
        tumor2count_inframe = {}
        total_count = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0}
        total_count_inframe = {"deletion": 0, "tandem_duplication": 0}
        for tumor_type in tumor_type_list.keys(): tumor2count[tumor_type] = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0} 
        for tumor_type in tumor_type_list.keys(): tumor2count_inframe[tumor_type] = {"deletion": 0, "tandem_duplication": 0}
        for tumor_type, sample, var, is_inframe in tumor_sample_var:
            tumor2count[tumor_type][var] = tumor2count[tumor_type][var] + 1
            total_count[var] = total_count[var] + 1
            if is_inframe == True:
                tumor2count_inframe[tumor_type][var] = tumor2count_inframe[tumor_type][var] + 1
                total_count_inframe[var] = total_count_inframe[var] + 1

        count_bar = str(total_count["deletion"] + total_count["tandem_duplication"] + total_count["inversion"] + total_count["translocation"]) + '\t' + \
                        str(total_count["deletion"]) + ',' + str(total_count["tandem_duplication"]) + ',' + str(total_count["inversion"]) + ',' + str(total_count["translocation"])
        if args.inframe_info == True:
            count_bar = count_bar + '\t' + str(total_count_inframe["deletion"] + total_count_inframe["tandem_duplication"]) + '\t' + \
                            str(total_count_inframe["deletion"]) + ',' + str(total_count_inframe["tandem_duplication"])
        
        for tumor_type in sorted(tumor2count):
            count_bar = count_bar + '\t' + str(tumor2count[tumor_type]["deletion"] + tumor2count[tumor_type]["tandem_duplication"] + tumor2count[tumor_type]["inversion"] + tumor2count[tumor_type]["translocation"]) + '\t' + \
                        str(tumor2count[tumor_type]["deletion"]) + ',' + str(tumor2count[tumor_type]["tandem_duplication"]) + ',' + str(tumor2count[tumor_type]["inversion"]) + ',' + str(tumor2count[tumor_type]["translocation"])
            if args.inframe_info == True:
                count_bar = count_bar + '\t' +  str(tumor2count_inframe[tumor_type]["deletion"] + tumor2count_inframe[tumor_type]["tandem_duplication"]) + '\t' + \
                            str(tumor2count_inframe[tumor_type]["deletion"]) + ',' + str(tumor2count_inframe[tumor_type]["tandem_duplication"])

        info = gene2info[gene] if gene in gene2info else "---" + '\t' + "---" + '\t' + "---"

        print >> hout, gene + '\t' + count_bar + '\t' + info


    hout.close()


def filter_main(args):


    if not os.path.exists(args.result_file):
        raise ValueError("file not exists: " + args.result_file)

    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"
    
    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
    
    hout = open(args.output, 'w')

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    sv_good_list = utils.filter_sv_list(args.result_file, args.fisher_thres, args.tumor_freq_thres, args.normal_freq_thres,
                                         args.normal_depth_thres, args.inversion_size_thres, args.max_size_thres,
                                         args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_bed, ens_junc_bed, grch2ucsc, 
                                         control_tb, args.control_num_thres, False)

    if len(sv_good_list) > 0:
        for i in range(0, len(sv_good_list)):
            if args.closest_exon == True:
                chr_ucsc1 = grch2ucsc[sv_good_list[i][0]] if sv_good_list[i][0] in grch2ucsc else sv_good_list[i][0] 
                chr_ucsc2 = grch2ucsc[sv_good_list[i][3]] if sv_good_list[i][3] in grch2ucsc else sv_good_list[i][3]
                dist_to_exon, target_exon = utils.distance_to_closest(chr_ucsc1, sv_good_list[i][1], chr_ucsc2, sv_good_list[i][4], ref_exon_tb)
                if len(target_exon) == 0: target_exon = ["---"]
                print >> hout, '\t'.join(sv_good_list[i]) + '\t' + str(dist_to_exon) + '\t' + ';'.join(target_exon)
            else:
                print >> hout, '\t'.join(sv_good_list[i])

    hout.close()


def concentrate_main(args):

    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    hout = open(args.output + ".tmp.sv_list", 'w')

    tumor_type_list = {}
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres, args.normal_freq_thres,
                                                 args.normal_depth_thres, args.inversion_size_thres, args.max_size_thres,
                                                 args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_bed, ens_junc_bed, grch2ucsc, 
                                                 control_tb, args.control_num_thres, False)

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


def merge_control_main(args):

    annotation_dir = args.annotation_dir
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    # make directory for output if necessary
    if os.path.dirname(args.output_prefix) != "" and not os.path.exists(os.path.dirname(args.output_prefix)):
        os.makedirs(os.path.dirname(args.output_prefix))

    hout = open(args.output_prefix + ".tmp.bedpe", 'w')

    tumor_type_list = {}
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            sample, tumor_type, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            sv_good_list = utils.filter_sv_list(result_file, args.fisher_thres, args.tumor_freq_thres, args.normal_freq_thres,
                                                 args.normal_depth_thres, args.inversion_size_thres, args.max_size_thres,
                                                 args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_bed, ens_junc_bed, grch2ucsc, 
                                                 control_tb, args.control_num_thres, True)


            for i in range(0, len(sv_good_list)):
                print >> hout, sv_good_list[i][0] + '\t' + str(int(sv_good_list[i][1]) - 1) + '\t' + str(int(sv_good_list[i][1])) + '\t' + \
                                sv_good_list[i][3] + '\t' + str(int(sv_good_list[i][4]) - 1) + '\t' + str(int(sv_good_list[i][4])) + '\t' + \
                                sample + '\t' + sv_good_list[i][6] + '\t' + sv_good_list[i][2] + '\t' + sv_good_list[i][5]
 
    hout.close()

    hout = open(args.output_prefix + ".bedpe", 'w')
    subprocess.call(["sort", "-k1,1", "-k3,3n", "-k4,4", "-k6,6n", args.output_prefix + ".tmp.bedpe"], stdout = hout)
    hout.close()

    # compress and index
    subprocess.call(["bgzip", "-f", args.output_prefix + ".bedpe"])
    subprocess.call(["tabix", "-p", "bed", args.output_prefix + ".bedpe.gz"])

    # remove intermediate file
    subprocess.call(["rm", "-rf", args.output_prefix + ".tmp.bedpe"])

