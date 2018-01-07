#! /usr/bin/env python

import subprocess
import my_seq


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

        header_line = "#"
        while header_line.startswith("#"):
            header_line = hin.readline().rstrip('\n')

        header = header_line.split('\t')

        for i in range(0, len(header)):
            if header[i] == "Ref": ref_ind = i
            if header[i] == "Alt": alt_ind = i
            if header[i] == "readPairNum_tumor": tum_ref_ind = i
            if header[i] == "variantPairNum_tumor": tum_var_ind = i
            if header[i] == "readPairNum_normal": nor_ref_ind = i
            if header[i] == "variantPairNum_normal": nor_var_ind = i
            if header[i] == "P-value(fisher_realignment)": fisher_ind = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if len(F[ref_ind]) >= 10 or len(F[alt_ind]) >= 10:

                if F[nor_ref_ind] == "---": continue
                if F[nor_var_ind] == "---": continue
                if F[tum_ref_ind] == "---": continue
                if F[tum_var_ind] == "---": continue

                bed_key = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + F[2]
                read_info = F[tum_ref_ind] + '\t' + F[tum_var_ind] + '\t' + \
                                str(round(float(F[tum_var_ind]) / (float(F[tum_ref_ind]) + float(F[tum_var_ind])), 3)) + '\t' + \
                                F[nor_ref_ind] + '\t' + F[nor_var_ind] + '\t' + \
                                str(round(float(F[nor_var_ind]) / (float(F[nor_ref_ind]) + float(F[nor_var_ind])), 3)) + '\t' + F[fisher_ind]

                var_info = ""
                gene_annotation = "---" + '\t' + "---" + '\t' + "---" + '\t' + "---"
                other_info = "---" + '\t' + "0" + '\t' + "---" + '\t' + "---"
                # deletion
                if len(F[ref_ind]) >= 10:
                    var_info = F[0] + '\t' + str(int(F[1]) - 1) + '\t' + "+" + '\t' + \
                                   F[0] + '\t' + str(int(F[1]) + len(F[ref_ind])) + '\t' + "-" + '\t' + "---" + '\t' + "deletion"
                    print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info + '\t' + other_info


                # tandem_duplication
                elif len(F[alt_ind]) >= 10:
                    # tandem_duplication check
                    flanking_seq_1 = my_seq.get_seq(reference, F[0], int(F[1]) + 1, int(F[1]) + len(F[alt_ind]))
                    flanking_seq_1_match, coord1 = my_seq.exact_alignment(F[alt_ind], flanking_seq_1)
                    flanking_seq_2 = my_seq.get_seq(reference, F[0], int(F[1]) - len(F[alt_ind]) + 1, int(F[1])) 
                    flanking_seq_2_match, coord2 = my_seq.exact_alignment(F[alt_ind], flanking_seq_2)
    
                    if flanking_seq_1_match == len(F[alt_ind]) or flanking_seq_2_match == len(F[alt_ind]):
                        var_info = F[0] + '\t' +  str(int(F[1]) + 1) + '\t' + "-" + '\t' + \
                                       F[0] + '\t' + str(int(F[1]) + len(F[alt_ind])) + '\t' + "+" + '\t' + "---" + '\t' + "tandem_duplication"
                        print >> hout, bed_key + '\t' + var_info + '\t' + gene_annotation + '\t' + read_info + '\t' + other_info
 
    hout.close()

    subprocess.call(["bgzip", "-f", output_file_prefix + ".bed"])
    subprocess.call(["tabix", "-p", "bed", output_file_prefix + ".bed.gz"])


