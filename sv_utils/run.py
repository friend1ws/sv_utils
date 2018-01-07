#! /usr/bin/env python

import sys, os, re, subprocess, gzip
import pysam
import utils, my_seq
from header_info import *
from utils import header_check

def count_main(args):

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

            print >> sys.stderr, "reading: " + result_file

            if args.inseq == True:
                type2count = {"deletion_nonseq": 0, "deletion_inseq": 0, "tandem_duplication_nonseq": 0, "tandem_duplication_inseq": 0,
                                "inversion_nonseq": 0, "inversion_inseq": 0, "translocation_nonseq": 0, "translocation_inseq": 0}

                with open(result_file, 'r') as hin:
                    for line in hin:
               
                        if line.startswith("#"): continue 
                        if utils.header_check(line.rstrip('\n')):
                            line = line.rstrip('\n')
                            header_info.read(line)
                            continue

                        F = line.rstrip('\n').split('\t')
                        if F[header_info.inserted_seq] != "---":
                            type2count[F[header_info.variant_type] + "_inseq"] = type2count[F[header_info.variant_type] + "_inseq"] + 1
                        else:
                            type2count[F[header_info.variant_type] + "_nonseq"] = type2count[F[header_info.variant_type] + "_nonseq"] + 1


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

                with open(result_file, 'r') as hin:
                    for line in hin:

                        if line.startswith("#"): continue
                        if line.startswith("Chr_1" + '\t' + "Pos_1"):
                            line = line.rstrip('\n')
                            header_info.read(line)
                            continue

                        F = line.rstrip('\n').split('\t')
                        type2count[F[header_info.variant_type]] = type2count[F[header_info.variant_type]] + 1

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

            with open(result_file, 'r') as hin:
                for line in hin:

                    if line.startswith('#'): continue
                    if utils.header_check(line.rstrip('\n')):
                        header_info.read(line.rstrip('\n'))
                        continue

                    F = line.rstrip('\n').split('\t')
                    var_type = F[header_info.variant_type]
                    genes1 = F[header_info.gene_1].split(';') if F[header_info.gene_1] != "---" else []
                    genes2 = F[header_info.gene_2].split(';') if F[header_info.gene_2] != "---" else []

                    # check in-frame or not
                    is_inframe = False
                    if var_type == "deletion":
                        var_size = int(F[header_info.pos_2]) - int(F[header_info.pos_1]) - 1
                        var_size = var_size - (0 if F[header_info.inserted_seq] == "---" else len(F[header_info.inserted_seq]))
                        if var_size % 3 == 0: is_inframe = True
                    elif var_type == "tandem_duplication":
                        var_size = int(F[header_info.pos_2]) - int(F[header_info.pos_1]) + 1
                        var_size = var_size + (0 if F[header_info.inserted_seq] == "---" else len(F[header_info.inserted_seq]))
                        if var_size % 3 == 0: is_inframe = True

                    for gene in list(set(genes1 + genes2)):
                        if gene in gene2type_sample: 
                            gene2type_sample[gene] = gene2type_sample[gene] + [(tumor_type, sample, var_type, is_inframe)]
                        else:
                            gene2type_sample[gene] = [(tumor_type, sample, var_type, is_inframe)]


    if args.inframe_info == True:
        print >> hout, '\t'.join(["gene", "all_type_total", "all_type_each", "all_type_total_inframe", "all_type_each_inframe"] + \
                    reduce(lambda x, y: x + y, [[x + "_total", x + "_each", x + "_total_inframe", x + "_each_inframe"] for x in sorted(tumor_type_list.keys())]) + \
                    ["Lawrence et al", "CGC", "kinbase", "Ye et al"])
    else:
        print >> hout, '\t'.join(["gene", "all_type_total", "all_type_each"] + \
                    reduce(lambda x, y: x + y, [[x + "_total", x + "_each"] for x in sorted(tumor_type_list.keys())]) + ["Lawrence et al", "CGC", "kinbase", "Ye et al"])

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

        info = gene2info[gene] if gene in gene2info else "---" + '\t' + "---" + '\t' + "---" + '\t' + "---"

        print >> hout, gene + '\t' + count_bar + '\t' + info


    hout.close()


def filter_main(args):

    import filter

    if not os.path.exists(args.input_file):
        raise ValueError("file not exists: " + args.input_file)

    # make directory for output if necessary
    if os.path.dirname(args.output_file) != "" and not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))
    
    hout = open(args.output_file, 'w')

    sv_good_list = filter.filter_sv_list(args)

    
    for i in range(0, len(sv_good_list)):

        # for meta info print
        if sv_good_list[i][0].startswith("#"):
            print >> hout, sv_good_list[i]
            continue
 
        # for header print
        if sv_good_list[i][header_info.chr_1] == "Chr_1" and sv_good_list[i][header_info.pos_1] == "Pos_1":
            print_header = '\t'.join(sv_good_list[i])
            print >> hout, print_header
            continue

        print_line = '\t'.join(sv_good_list[i])

        print >> hout, print_line



def mutation_main(args):

    import mutation

    mutation.make_mut_db(args.mutation_result_file, args.output_file + ".mutation", args.reference)
    mut_tb = pysam.TabixFile(args.output_file + ".mutation.bed.gz")

    hout = open(args.output_file, 'w')

    dup_list = {}
    with open(args.sv_result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"):
                print >> hout, line.rstrip('\n')
                continue

            if header_check(line.rstrip('\n')):
                header_info.read(line.rstrip('\n'))
                print >> hout, line.rstrip('\n') + '\t' + "Mutation_Detection"
                continue

            F = line.rstrip('\n').split('\t')

            duplicated_flag = 0

            if F[header_info.variant_type] in ["deletion", "tandem_duplication"] and \
                abs(int(F[header_info.pos_1]) - int(F[header_info.pos_2])) <= 100:

                # check exon annotation for the side 1
                tabixErrorFlag = 0
                try:
                    records = mut_tb.fetch(F[header_info.chr_1], int(F[header_info.pos_1]) - 30, int(F[header_info.pos_2]) + 30)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorFlag = 1

                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        mut_length = int(record[7]) - int(record[4]) - 1
                        sv_length = int(F[header_info.pos_2]) - int(F[header_info.pos_1]) - 1

                        if float(abs(mut_length - sv_length)) / sv_length < 0.2 and record[10] == F[header_info.variant_type]:

                            duplicated_flag = 1
                            dup_list[record[0] + '\t' + record[1] + '\t' + record[2]] = 1

            if duplicated_flag == 1:
                print >> hout, '\t'.join(F) + '\t' + "mut;sv"
            else:
                print >> hout, '\t'.join(F) + '\t' + "sv"


    with gzip.open(args.output_file + ".mutation.bed.gz") as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            bed_key = F[0] + '\t' + F[1] + '\t' + F[2]
            if bed_key in dup_list: continue
            print >> hout, '\t'.join(F) + "mut"

    subprocess.call(["rm", "-rf", args.output_file + ".mutation.bed.gz"])
    subprocess.call(["rm", "-rf", args.output_file + ".mutation.bed.gz.tbi"])

    hout.close()


def concentrate_main(args):

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

            with open(result_file, 'r') as hin:
                for line in hin:
                    if utils.header_check(line.rstrip('\n')):
                        line = line.rstrip('\n')
                        header_info.read(line)
                        continue

                    F = line.rstrip('\n').split('\t')
                    if F[7] in ["deletion", "tandem_duplication"]:
                        print >> hout, sample + '\t' + tumor_type + '\t' + '\t'.join(F)

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

    import genomonsv.mergeFunction, genomonsv.utils

    # make directory for output if necessary
    if os.path.dirname(args.output_file) != "" and not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    hout = open(args.output_file + ".temp", 'w')

    tumor_type_list = {}
    gene2type_sample = {}
    with open(args.result_list, 'r') as hin:
        for line in hin:

            label, tumor_type, result_file = line.rstrip('\n').split('\t')
            # label, result_file = line.rstrip('\n').split('\t')
            if tumor_type not in tumor_type_list: tumor_type_list[tumor_type] = 1

            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)


            num = 1
            with open(result_file, 'r') as hin:
                for line in hin:
                    if line.startswith("#"): continue
                    if utils.header_check(line.rstrip('\n')):
                        line = line.rstrip('\n')
                        header_info.read(line)
                        continue

                    F = line.rstrip('\n').split('\t')
                    inseqLen = len(F[header_info.inserted_seq]) if F[header_info.inserted_seq] != "---" else 0

                    print >> hout, '\t'.join([F[header_info.chr_1], str(int(F[header_info.pos_1]) - 1), F[header_info.pos_1], \
                                              F[header_info.chr_2], str(int(F[header_info.pos_2]) - 1), F[header_info.pos_2], \
                                              "junction_" + str(num),  str(inseqLen), \
                                              F[header_info.dir_1], F[header_info.dir_2], label, "1"])

                    num = num + 1

    hout.close()

    """
    hout = open(args.output_prefix + ".bedpe", 'w')
    subprocess.call(["sort", "-k1,1", "-k3,3n", "-k4,4", "-k6,6n", args.output_prefix + ".tmp.bedpe"], stdout = hout)
    hout.close()

    # compress and index
    subprocess.call(["bgzip", "-f", args.output_prefix + ".bedpe"])
    subprocess.call(["tabix", "-p", "bed", args.output_prefix + ".bedpe.gz"])

    # remove intermediate file
    subprocess.call(["rm", "-rf", args.output_prefix + ".tmp.bedpe"])
    """

    # utils.processingMessage("sorting the aggregated junction file")
    genomonsv.utils.sortBedpe(args.output_file + ".temp", args.output_file + ".temp.sort")

    # utils.processingMessage("merging the same junction in the aggregated junction file")
    genomonsv.mergeFunction.organizeControl(args.output_file + ".temp.sort", args.output_file + ".temp.merged", 20)

    # utils.processingMessage("sorting the merged junction file")
    genomonsv.utils.sortBedpe(args.output_file + ".temp.merged", args.output_file + ".temp.merged.sort")

    # utils.processingMessage("compressing the merged junction file")
    genomonsv.utils.compress_index_bed(args.output_file + ".temp.merged.sort", args.output_file)



def realign_main(args):

    from genomonsv import filterFunction

    if args.tumor_bam is None:
        print >> sys.stderr, "tumor_bam file should be input"
        sys.exit(1)

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    matchedControlFlag = True if args.control_bam is not None else False

    # generate bedpe file
    hout = open(args.output + ".tmp1.bedpe", 'w')
    i = 0
    with open(args.result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                line = line.rstrip('\n')
                header_info.read(line)
                continue

            F = line.rstrip('\n').split('\t')
            print >> hout, '\t'.join([F[header_info.chr_1], str(int(F[header_info.pos_1]) - 1), F[header_info.pos_1], \
                                      F[header_info.chr_2], str(int(F[header_info.pos_2]) - 1), F[header_info.pos_2], \
                                      "genoemonSV_" + str(i), F[header_info.inserted_seq], F[header_info.dir_1], F[header_info.dir_2]] + \
                                      ["---" for i in range(14)])
            i = i + 1

    hout.close()


    filterFunction.validateByRealignment(args.output + ".tmp1.bedpe",
                    args.output + ".tmp2.bedpe",
                    args.tumor_bam,
                    args.control_bam,
                    args.reference,
                    "-stepSize=5 -repMatch=2253",
                    500,
                    5000,
                    1000,
                    5,
                    1000,
                    1000)

    key2AF_info = {}
    with open(args.output + ".tmp2.bedpe", 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            key = '\t'.join(F[:7])

            tumorAF = 0 
            if float(F[7]) + float(F[8]) > 0: tumorAF = float(F[8]) / (float(F[7]) + float(F[8]))     
            tumorAF = str(round(tumorAF, 4))

            normalAF = "---"
            if matchedControlFlag == True:
                normalAF = 0
                if float(F[9]) + float(F[10]) > 0: normalAF = float(F[10]) / (float(F[9]) + float(F[10]))
                normalAF = str(round(normalAF, 4))

            if matchedControlFlag == True:
                key2AF_info[key] = '\t'.join([F[7], F[8], tumorAF, F[9], F[10], normalAF, F[11]])
            else:
                key2AF_info[key] = '\t'.join([F[7], F[8], tumorAF])


    hout = open(args.output, 'w') 
    with open(args.result_file, 'r') as hin:
        for line in hin:
            
            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                line = line.rstrip('\n')
                if matchedControlFlag == True:
                    print >> hout, line + '\t' + "Num_Tumor_Ref_Read_Pair_re" + '\t' + "Num_Tumor_Var_Read_Pair_re" + '\t' + "Tumor_VAF_re" + '\t' + \
                                                 "Num_Control_Ref_Read_Pair_re" + '\t'+ "Num_Control_Var_Read_Pair_re" + '\t' + "Control_VAF_re" + '\t' + \
                                                 "Minus_Log_Fisher_P_value_re" 
                else:
                    print >> hout, line + '\t' + "Num_Tumor_Ref_Read_Pair_re" + '\t' + "Num_Tumor_Var_Read_Pair_re" + '\t' + "Tumor_VAF_re"
                continue

            F = line.rstrip('\n').split('\t')
            key = '\t'.join(F[:7])
            print >> hout, '\t'.join(F) + '\t' + key2AF_info[key]

    hout.close()

    subprocess.call(["rm", "-rf", args.output + ".tmp1.bedpe"])
    subprocess.call(["rm", "-rf", args.output + ".tmp2.bedpe"])


def primer_main(args):
 
    from genomonsv import realignmentFunction
    from primer3 import bindings 

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    param = {"reference_genome": args.reference, "split_refernece_thres": 1000, "validate_sequence_length": 250}

    possible_chr = [str(x) for x in range(1, 23)] + ['X', 'Y'] + \
                   ["chr" + str(x) for x in range(1, 23)] + ['chrX', 'chrY']

    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:
    
            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                line = line.rstrip('\n')
                header_info.read(line)
                print >> hout, line + '\t' + "Primer1" + '\t' + "Primer2" + '\t' + "Primer3" + '\t' + "Primer4" + '\t' + "Primer5"
                continue
            
            F = line.rstrip('\n').split('\t')
            chr1, pos1, dir1, chr2, pos2, dir2, junc_seq = F[header_info.chr_1], F[header_info.pos_1], F[header_info.dir_1], \
                                                           F[header_info.chr_2], F[header_info.pos_2], F[header_info.dir_2], F[header_info.inserted_seq]

            if chr1 not in possible_chr or chr2 not in possible_chr:
                print >> sys.stderr, "Skip a SV incolving atypical chromosomes: %s,%s,%s,%s,%s,%s" % \
                            (chr1, pos1, dir1, chr2, pos2, dir2)
                continue
            
            junc_seq_len = 0 if junc_seq == "---" else len(junc_seq)

            realignmentFunction.getRefAltForSV(args.output + ".contig.tmp.fa", chr1, pos1, dir1, chr2, pos2, dir2, junc_seq, args.reference, 1000, 250)

            with open(args.output + ".contig.tmp.fa") as hin2:
                lines2 = hin2.readlines()
                for i in range(len(lines2)):
                    lines2[i] = lines2[i].rstrip('\n')
                    if lines2[i].startswith('>') and lines2[i].endswith("alt"):
                        seq = lines2[i + 1].rstrip('\n')

                        primer = bindings.designPrimers(
                            {
                                'SEQUENCE_ID': 'MH1000',
                                'SEQUENCE_TEMPLATE': seq,
                                'SEQUENCE_TARGET': [225,50 + junc_seq_len],
                                'SEQUENCE_INCLUDED_REGION': [10, len(seq) - 20]
                            },
                            {
                                'PRIMER_PRODUCT_SIZE_RANGE': [[150,250],[100,300],[301,400],[401,500]],
                            })

                        primer_left_right = ["---"] * 5
                        for i in range(5):
                            if "PRIMER_LEFT_" + str(i) + "_SEQUENCE" in primer and "PRIMER_RIGHT_" + str(i) + "_SEQUENCE" in primer and \
                               "PRIMER_LEFT_" + str(i) + "_TM" in primer and "PRIMER_RIGHT_" + str(i) + "_TM" in primer and \
                               "PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE" in primer:
                                primer_left_right[i] = primer["PRIMER_LEFT_" + str(i) + "_SEQUENCE"] + ";" + primer["PRIMER_RIGHT_" + str(i) + "_SEQUENCE"] + ';' + \
                                                       str(round(primer["PRIMER_LEFT_" + str(i) + "_TM"], 3)) + ";" + str(round(primer["PRIMER_RIGHT_" + str(i) + "_TM"], 3)) + ';' + \
                                                       str(primer["PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"])
 
                        print >> hout, '\t'.join(F) + '\t' + '\t'.join(primer_left_right)
              

    hout.close()    
    subprocess.call(["rm", "-rf", args.output + ".contig.tmp.fa"])
 

def vcf_main(args):

    # generate bedpe file
    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:
    
            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                header_info.read(line.rstrip('\n'))
                continue

            F = line.rstrip('\n').split('\t')

            if F[header_info.variant_type] in ["inversion", "translocation"]: continue
            if abs(int(F[header_info.pos_1]) - int(F[header_info.pos_2])) > int(args.max_size_thres): continue

            if F[header_info.variant_type] == "deletion":
                ref_seq = my_seq.get_seq(args.reference, F[header_info.chr_1], int(F[header_info.pos_1]), int(F[header_info.pos_2]) - 1)
                alt_seq = ref_seq[0] if F[header_info.inserted_seq] == "---" else ref_seq[0] + F[header_info.inserted_seq] 
                pos = F[header_info.pos_1]
            elif F[header_info.variant_type] == "tandem_duplication":
                alt_seq = my_seq.get_seq(args.reference, F[header_info.chr_1], int(F[header_info.pos_1]) - 1, int(F[header_info.pos_2]))
                alt_seq = alt_seq if F[header_info.inserted_seq] == "---" else alt_seq + F[header_info.inserted_seq] 
                ref_seq = alt_seq[0]
                pos = str(int(F[header_info.pos_1]) - 1)

            print >> hout, '\t'.join([F[header_info.chr_1], pos, '.', ref_seq, alt_seq, '.', "PASS", '.']) 

    hout.close()


def homology_main(args):

    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                header_info.read(line.rstrip('\n'))
                print_header = line.rstrip('\n') + '\t' + "Homology_Match"
                print >> hout, print_header
                continue

            F = line.rstrip('\n').split('\t')

            # for meta info print
            if F[0].startswith("#"):
                print >> hout, '\t'.join(F)
                continue

            var_size = 500000
            if F[header_info.variant_type] == "deletion":
                var_size = int(F[header_info.pos_2]) - int(F[header_info.pos_1]) - 1
            elif F[header_info.variant_type] == "tandem_duplication":
                var_size = int(F[header_info.pos_2]) - int(F[header_info.pos_1]) + 1
            
            homology_match = utils.check_homology(F[header_info.chr_1], F[header_info.pos_1], F[header_info.dir_1],
                                                  F[header_info.chr_2], F[header_info.pos_2], F[header_info.dir_2], 
                                                  args.reference, min(var_size, 100))

            print >> hout, '\t'.join(F) + '\t' + str(homology_match)
 
    hout.close()


def nonB_DB_main(args):

    all_nonB_DB_type = ["A_Phased_Repeat", "Direct_Repeat", "G_Quadruplex_Motif", "Inverted_Repeat", 
                       "Mirror_Repeat", "Short_Tandem_Repeat", "Z_DNA_Motif"]

    if not os.path.exists(args.result_file):
        raise ValueError("file not exists: " + args.result_file)

    annotation_dir = args.annotation_dir
    nonB_DB_bed = annotation_dir + "/nonB_DB.bed.gz"
    nonB_DB_tb = pysam.TabixFile(nonB_DB_bed)

    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]


    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:
            if utils.header_check(line.rstrip('\n')):
                header_info.read(line.rstrip('\n'))
                print_header = line.rstrip('\n') + '\t' + '\t'.join([x + "_dist1" + '\t' + x + "_dist2" for x in all_nonB_DB_type])
                print >> hout, print_header
                continue
            
            F = line.rstrip('\n').split('\t')

            # for meta info print
            if F[0].startswith("#"):
                print >> hout, '\t'.join(F)
                continue

            chr_ucsc1 = grch2ucsc[F[header_info.chr_1]] if F[header_info.chr_1] in grch2ucsc else F[header_info.chr_1]
            chr_ucsc2 = grch2ucsc[F[header_info.chr_2]] if F[header_info.chr_2] in grch2ucsc else F[header_info.chr_2]

            print_dist_bar = ''
            for nonB_DB_type in all_nonB_DB_type:
                nonB_DB_dist1 = utils.nonB_DB_dist_check(chr_ucsc1, int(F[header_info.pos_1]), F[header_info.dir_1], nonB_DB_tb, nonB_DB_type)
                nonB_DB_dist2 = utils.nonB_DB_dist_check(chr_ucsc2, int(F[header_info.pos_2]), F[header_info.dir_2], nonB_DB_tb, nonB_DB_type)
                print_dist_bar = print_dist_bar + '\t' + str(nonB_DB_dist1) + '\t' + str(nonB_DB_dist2)

            print >> hout, '\t'.join(F) + print_dist_bar



def RSS_main(args):

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    rss_pwm = my_seq.generate_rss_pwm()

    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                line = line.rstrip('\n')
                header_info.read(line)
                print >> hout, line + '\t' + "RSS_score_1" + '\t' + "RSS_info_1" + '\t' + "RSS_score_2" + '\t' + "RSS_info_2"
                continue

            F = line.rstrip('\n').split('\t')
            
            seq1 = my_seq.get_seq(args.reference, F[header_info.chr_1], int(F[header_info.pos_1]) - 50, int(F[header_info.pos_1]) + 50)
            seq2 = my_seq.get_seq(args.reference, F[header_info.chr_2], int(F[header_info.pos_2]) - 50, int(F[header_info.pos_2]) + 50)

            rss_info_1 = my_seq.get_max_rss_score(seq1, rss_pwm[0], rss_pwm[1])
            rss_info_2 = my_seq.get_max_rss_score(seq2, rss_pwm[0], rss_pwm[1])
   
            print >> hout, '\t'.join(F) + '\t' + \
                            str(round(rss_info_1[0], 3)) + '\t' + \
                            ';'.join([rss_info_1[1], rss_info_1[2], str(int(rss_info_1[3] - 50)), str(rss_info_1[4]), str(rss_info_1[5])]) + '\t' + \
                            str(round(rss_info_2[0], 3)) + '\t' + \
                            ';'.join([rss_info_2[1], rss_info_2[2], str(int(rss_info_2[3] - 50)), str(rss_info_2[4]), str(rss_info_2[5])])

    hout.close()



def AID_main(args):

    # make directory for output if necessary
    if os.path.dirname(args.output) != "" and not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    hout = open(args.output, 'w')
    with open(args.result_file, 'r') as hin:
        for line in hin:

            if line.startswith("#"): continue
            if utils.header_check(line.rstrip('\n')):
                line = line.rstrip('\n')
                header_info.read(line)
                print >> hout, line + '\t' + "CG_motif_info_1" + '\t' + "CG_motif_info_2" + '\t' + "WGCW_motif_info_1" + '\t' + "WGCW_motif_info_2"
                continue

            F = line.rstrip('\n').split('\t')

            seq1 = my_seq.get_seq(args.reference, F[header_info.chr_1], int(F[header_info.pos_1]) - args.check_size, int(F[header_info.pos_1]) + args.check_size)
            seq2 = my_seq.get_seq(args.reference, F[header_info.chr_2], int(F[header_info.pos_2]) - args.check_size, int(F[header_info.pos_2]) + args.check_size)

           
            CG_starts_1 = [match.start() - 10 for match in re.finditer(r'CG', seq1)]
            CG_starts_2 = [match.start() - 10 for match in re.finditer(r'CG', seq2)]
            WGCW_starts_1 = [match.start() - 10 for match in re.finditer(r'[AT]GC[AT]', seq1)]
            WGCW_starts_2 = [match.start() - 10 for match in re.finditer(r'[AT]GC[AT]', seq2)]

            if len(CG_starts_1) == 0: CG_starts_1.append("---")
            if len(CG_starts_2) == 0: CG_starts_2.append("---")
            if len(WGCW_starts_1) == 0: WGCW_starts_1.append("---")
            if len(WGCW_starts_2) == 0: WGCW_starts_2.append("---")

            print >> hout, '\t'.join(F) + '\t' + \
                            ','.join([str(x) for x in CG_starts_1]) + '\t' + \
                            ','.join([str(x) for x in CG_starts_2]) + '\t' + \
                            ','.join([str(x) for x in WGCW_starts_1]) + '\t' + \
                            ','.join([str(x) for x in WGCW_starts_2]) 

    hout.close()


