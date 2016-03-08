#! /usr/bin/env python

import sys, os, subprocess, gzip
import pysam
import utils

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

            if args.inseq == True:
                type2count = {"deletion_nonseq": 0, "deletion_inseq": 0, "tandem_duplication_nonseq": 0, "tandem_duplication_inseq": 0,
                                "inversion_nonseq": 0, "inversion_inseq": 0, "translocation_nonseq": 0, "translocation_inseq": 0}

                with open(result_file, 'r') as hin:
                    for line in hin:
                        F = line.rstrip('\n').split('\t')
                        if F[6] != "---":
                            type2count[F[7] + "_inseq"] = type2count[F[7] + "_inseq"] + 1
                        else:
                            type2count[F[7] + "_nonseq"] = type2count[F[7] + "_nonseq"] + 1


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
                        F = line.rstrip('\n').split('\t')
                        type2count[F[7]] = type2count[F[7]] + 1

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
                    F = line.rstrip('\n').split('\t')
                    var_type = F[7]
                    genes1 = F[8].split(';') if F[8] != "---" else []
                    genes2 = F[9].split(';') if F[9] != "---" else []

                    # check in-frame or not
                    is_inframe = False
                    if var_type == "deletion":
                        var_size = int(F[4]) - int(F[1]) - 1
                        var_size = var_size - (0 if F[6] == "---" else len(F[6]))
                        if var_size % 3 == 0: is_inframe = True
                    elif var_type == "tandem_duplication":
                        var_size = int(F[4]) - int(F[1]) + 1
                        var_size = var_size + (0 if F[6] == "---" else len(F[6]))
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

    if args.mutation_result != "" and args.reference == "":
        print >> sys.stderr, "input reference when using mutation_result"
        sys.exit(1)

    annotation_dir = args.annotation_dir
    ref_gene_bed = annotation_dir + "/refGene.bed.gz"
    ens_gene_bed = annotation_dir + "/ensGene.bed.gz"
    ref_junc_bed = annotation_dir + "/refJunc.bed.gz"
    ens_junc_bed = annotation_dir + "/ensJunc.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    ref_coding_bed = annotation_dir + "/refCoding.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"
    simple_repeat_bed = annotation_dir + "/simpleRepeat.bed.gz"
 
    ref_gene_tb = pysam.TabixFile(ref_gene_bed)
    ens_gene_tb = pysam.TabixFile(ens_gene_bed)
    ref_junc_tb = pysam.TabixFile(ref_junc_bed)
    ens_junc_tb = pysam.TabixFile(ens_junc_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)
    ref_coding_tb = pysam.TabixFile(ref_coding_bed)
    control_tb = pysam.TabixFile(args.control) if args.control is not None else None
    simple_repeat_tb = pysam.TabixFile(simple_repeat_bed) if args.remove_simple_repeat is not None else None

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
                                        args.within_exon, ref_exon_tb, ens_exon_tb, ref_junc_tb, ens_junc_tb,
                                        simple_repeat_tb, grch2ucsc, control_tb, args.control_num_thres, False)

    mut_tb = None
    if args.mutation_result != "":
        utils.make_mut_db(args.mutation_result, args.output + ".mutation", args.reference)
        mut_tb = pysam.TabixFile(args.output + ".mutation.bed.gz")

    dup_list = {}
    for i in range(0, len(sv_good_list)):
        chr_ucsc1 = grch2ucsc[sv_good_list[i][0]] if sv_good_list[i][0] in grch2ucsc else sv_good_list[i][0]
        chr_ucsc2 = grch2ucsc[sv_good_list[i][3]] if sv_good_list[i][3] in grch2ucsc else sv_good_list[i][3]

        if args.re_annotation == True:
            sv_good_list[i][8], sv_good_list[i][9], sv_good_list[i][10], sv_good_list[i][11] = \
                utils.get_gene_annotation(chr_ucsc1, sv_good_list[i][1], chr_ucsc2, sv_good_list[i][4], ref_gene_tb, ref_exon_tb)

        print_line = '\t'.join(sv_good_list[i])

        if args.mutation_result != "":
            if sv_good_list[i][7] in ["deletion", "tandem_duplication"] and abs(int(sv_good_list[i][1]) - int(sv_good_list[i][4])) <= 100:

                # check exon annotation for the side 1
                tabixErrorFlag = 0
                try:
                    records = mut_tb.fetch(sv_good_list[i][0], int(sv_good_list[i][1]) - 50, int(sv_good_list[i][4]) + 50)
                except Exception as inst:
                    # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                    tabixErrorFlag = 1

                duplicated_flag = 0
                if tabixErrorFlag == 0:
                    for record_line in records:
                        record = record_line.split('\t')
                        if int(record[7]) - int(record[4]) == int(sv_good_list[i][4]) - int(sv_good_list[i][1]) and record[10] == sv_good_list[i][7]:
                            duplicated_flag = 1
                            dup_list[record[0] + '\t' + record[1] + '\t' + record[2]] = 1

                print_line = print_line + '\t' + "mut;sv" if duplicated_flag == 1 else print_line + '\t' + "sv"
            else:
                print_line = print_line + '\t' + "sv"   
            
        if args.closest_exon == True:
            dist_to_exon, target_exon = utils.distance_to_closest(chr_ucsc1, sv_good_list[i][1], chr_ucsc2, sv_good_list[i][4], ref_exon_tb)
            if len(target_exon) == 0: target_exon = ["---"]
            print_line = print_line  + '\t' + str(dist_to_exon) + '\t' + ';'.join(target_exon)

        if args.closest_coding == True:
            dist_to_exon, target_exon = utils.distance_to_closest(chr_ucsc1, sv_good_list[i][1], chr_ucsc2, sv_good_list[i][4], ref_coding_tb, False)
            if len(target_exon) == 0: target_exon = ["---"]
            print_line = print_line  + '\t' + str(dist_to_exon) + '\t' + ';'.join(target_exon)
 
        if args.coding_info == True:
            # within gene or accross gene ?
            gene_flag = 0
            within_gene_flag = 0
            for (g1, g2) in zip(sv_good_list[i][8].split(';'), sv_good_list[i][9].split(';')):
                if g1 != "---" and g1 == g2: within_gene_flag = 1
                if g1 != "---" or g2 != "---": gene_flag = 1

            if within_gene_flag == 1:
                coding_info = utils.check_coding_info(chr_ucsc1, sv_good_list[i][1], sv_good_list[i][4], ref_coding_tb)
                print_line = print_line + '\t' + "within_gene" + '\t' + coding_info                
            elif gene_flag == 1:
                print_line = print_line + '\t' + "across_gene" + '\t' + "---\t---" 
            else:
                print_line = print_line + '\t' + "intergenic" + '\t' + "---\t---"

        if args.fusion_info is not None:
            print_line = print_line + '\t' + utils.check_fusion_direction(chr_ucsc1, sv_good_list[i][1], sv_good_list[i][2],
                                                                          chr_ucsc2, sv_good_list[i][4], sv_good_list[i][5], ref_gene_tb)

        print >> hout, print_line


    if args.mutation_result != "":
        with gzip.open(args.output + ".mutation.bed.gz") as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                bed_key = F[0] + '\t' + F[1] + '\t' + F[2]
                if bed_key in dup_list: continue

                chr_ucsc1 = grch2ucsc[F[3]] if F[3] in grch2ucsc else F[3] 
                chr_ucsc2 = grch2ucsc[F[6]] if F[6] in grch2ucsc else F[6] 
                F[11], F[12], F[13], F[14] = utils.get_gene_annotation(chr_ucsc1, F[4], chr_ucsc2, F[7], ref_gene_tb, ref_exon_tb)
                print_line = '\t'.join(F[3:]) + '\t' + "mut"

                if args.closest_exon == True:
                    dist_to_exon, target_exon = utils.distance_to_closest(chr_ucsc1, F[4], chr_ucsc2, F[7], ref_exon_tb)
                    if len(target_exon) == 0: target_exon = ["---"]
                    print_line = print_line  + '\t' + str(dist_to_exon) + '\t' + ';'.join(target_exon)

                if args.closest_coding == True:
                    dist_to_exon, target_exon = utils.distance_to_closest(chr_ucsc1, F[4], chr_ucsc2, F[7], ref_coding_tb, False)
                    if len(target_exon) == 0: target_exon = ["---"]
                    print_line = print_line  + '\t' + str(dist_to_exon) + '\t' + ';'.join(target_exon)

                if args.coding_info == True:
                    # within gene or accross gene ?
                    gene_flag = 0
                    within_gene_flag = 0
                    for (g1, g2) in zip(F[11].split(';'), F[12].split(';')):
                        if g1 != "---" and g1 == g2: within_gene_flag = 1
                        if g1 != "---" or g2 != "---": gene_flag = 1

                    if within_gene_flag == 1:
                        chr_ucsc1 = grch2ucsc[F[3]] if F[3] in grch2ucsc else F[3]
                        coding_info = utils.check_coding_info(chr_ucsc1, F[4], F[7], ref_coding_tb)
                        print_line = print_line + '\t' + "within_gene" + '\t' + coding_info
                    elif gene_flag == 1:
                        print_line = print_line + '\t' + "across_gene" + '\t' + "---\t---"
                    else:
                        print_line = print_line + '\t' + "intergenic" + '\t' + "---\t---"

                if args.fusion_info is not None:
                    print_line = print_line + '\t' + "---"

                print >> hout, print_line

        subprocess.call(["rm", "-rf", args.output + ".mutation.bed.gz"])
        subprocess.call(["rm", "-rf", args.output + ".mutation.bed.gz.tbi"])

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

            with open(result_file, 'r') as hin:
                for line in hin:
                    F = line.rstrip('\n').split('\t')
                    print >> hout, F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[1] + '\t' + F[3] + '\t' + str(int(F[4]) - 1) + '\t' + F[4] + '\t' + \
                                    sample + '\t' + F[6] + '\t' + F[2] + '\t' + F[5]
    hout.close()

    hout = open(args.output_prefix + ".bedpe", 'w')
    subprocess.call(["sort", "-k1,1", "-k3,3n", "-k4,4", "-k6,6n", args.output_prefix + ".tmp.bedpe"], stdout = hout)
    hout.close()

    # compress and index
    subprocess.call(["bgzip", "-f", args.output_prefix + ".bedpe"])
    subprocess.call(["tabix", "-p", "bed", args.output_prefix + ".bedpe.gz"])

    # remove intermediate file
    subprocess.call(["rm", "-rf", args.output_prefix + ".tmp.bedpe"])



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
            F = line.rstrip('\n').split('\t')
            print >> hout, '\t'.join([F[0], str(int(F[1]) - 1), F[1], F[3], str(int(F[4]) - 1), F[4], "genoemonSV_" + str(i), F[6], F[2], F[5]])
            i = i + 1

    hout.close()

    # yaml input
    param = {"max_depth": 5000, "search_length": 1000, "search_margin": 5, "reference_genome": args.reference,
             "split_refernece_thres": 1000, "validate_sequence_length": 1000, "STD_thres": 500}

    filterFunction.validateByRealignment(args.output + ".tmp1.bedpe",
                    args.output + ".tmp2.bedpe",
                    args.tumor_bam,
                    args.control_bam,
                    "blat -stepSize=5 -repMatch=2253",
                    matchedControlFlag,
                    param)

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
            F = line.rstrip('\n').split('\t')
            key = '\t'.join(F[:7])
            print >> hout, '\t'.join(F) + '\t' + key2AF_info[key]

    hout.close()

    subprocess.call(["rm", "-rf", args.output + ".tmp1.bedpe"])
    subprocess.call(["rm", "-rf", args.output + ".tmp2.bedpe"])

