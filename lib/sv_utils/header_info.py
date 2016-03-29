#! /usr/bin/env python

class Header_info(object):

    def __init__(self):

        self.chr_1 = 0
        self.pos_1 = 1
        self.dir_1 = 2
        self.chr_2 = 3
        self.pos_2 = 4
        self.dir_2 = 5
        self.inserted_seq = 6
        self.variant_type = 7
        self.gene_1 = 8
        self.gene_2 = 9
        self.exon_1 = 10
        self.exon_2 = 11
        self.num_tumor_ref_read_pair = 12
        self.num_tumor_var_read_pair = 13
        self.tumor_vaf = 14
        self.num_control_ref_read_pair = 15
        self.num_control_var_read_pair = 16
        self.control_vaf = 17
        self.minus_log_fisher_p_value = 18
        self.non_matched_control_sample_with_max_junction = 19
        self.num_non_matched_control_junction = 20
        self.max_over_hang_1 = 21
        self.max_over_hang_2 = 22


    def read(self, header):
    
        F = header.rstrip('\n').split('\t')
        for i in range(len(F)):
            if F[i] == "Chr_1": self.chr_1 = i
            if F[i] == "Pos_1": self.pos_1 = i
            if F[i] == "Dir_1": self.dir_1 = i
            if F[i] == "Chr_2": self.chr_2 = i
            if F[i] == "Pos_2": self.pos_2 = i
            if F[i] == "Dir_2": self.dir_2 = i
            if F[i] == "Inserted_Seq": self.inserted_seq = i
            if F[i] == "Variant_Type": self.variant_type = i
            if F[i] == "Gene_1": self.gene_1 = i
            if F[i] == "Gene_2": self.gene_2 = i
            if F[i] == "Exon_1": self.exon_1 = i
            if F[i] == "Exon_2": self.exon_2 = i
            if F[i] == "Num_Tumor_Ref_Read_Pair": self.num_tumor_ref_read_pair = i
            if F[i] == "Num_Tumor_Var_Read_Pair": self.num_tumor_var_read_pair = i
            if F[i] == "Tumor_VAF": self.tumor_vaf = i
            if F[i] == "Num_Control_Ref_Read_Pair": self.num_control_ref_read_pair = i
            if F[i] == "Num_Control_Var_Read_Pair": self.num_control_var_read_pair = i
            if F[i] == "Control_VAF": self.control_vaf = i
            if F[i] == "Minus_Log_Fisher_P_value": self.minus_log_fisher_p_value = i
            if F[i] == "Non-Matched_Control_Sample_With_Max_Junction": self.non_matched_control_sample_with_max_junction = i
            if F[i] == "Num_Max_Non-Matched_Control_Junction": self.num_non_matched_control_junction = i
            if F[i] == "Max_Over_Hang_1": self.max_over_hang_1 = i
            if F[i] == "Max_Over_Hang_2": self.max_over_hang_2 = i


global header_info
header_info = Header_info()

