#! /usr/bin/env python

from run import *
import argparse

def create_parser():

    parser = argparse.ArgumentParser(prog = "sv_utils")
    parser.add_argument("--version", action = "version", version = "sv_utils-0.5.0")

    subparsers = parser.add_subparsers()

    ##########
    # count_summary 
    count = subparsers.add_parser("count",
                                  help = "Summarize the frequency of each variant type (deletion, tandem_duplication, inversion, translocation) for each sample")

    count.add_argument("result_list", metavar = "result_list.txt", default = None, type = str,
                        help = "1st column: sample IDs, 2nd column: tumor type, 3rd column: genomon SV result path")

    count.add_argument("output", metavar = "output.txt", default = None, type = str,
                        help = "the path to the output file")

    count.add_argument("--inseq", default = False, action = "store_true",
                       help = "separete variant by existence of inserted sequences")

    count.set_defaults(func = count_main)

    ##########
    # gene_summary_summary 
    gene_summary = subparsers.add_parser("gene_summary",
                                         help = "Summarize the frequency of each variant type (deletion, tandem_duplication, inversion, translocation) for each cancer gene")

    gene_summary.add_argument("result_list", metavar = "result_list.txt", default = None, type = str,
                              help = "1st column: sample IDs, 2nd column: tumor type, 3rd column: genomon SV result path")

    gene_summary.add_argument("output", metavar = "output.txt", default = None, type = str,
                              help = "the path to the output file")

    gene_summary.add_argument("--cancer_gene_list", metavar = "cancer_gene_list", default = None, type = str, required=True,
                              help = "the path to the cancer list file")

    gene_summary.add_argument("--inframe_info", default = False, action = "store_true",
                             help = "add inframe information")

    gene_summary.set_defaults(func = gene_summary_main)


    ##########
    # filter
    filter_parser = subparsers.add_parser("filter",
                                        help = "Filter out GenomonSV results outside specified conditions")

    filter_parser.add_argument("input_file", metavar = "genomonSV.result.txt", default = None, type = str,
                        help = "the path to genomon SV result")

    filter_parser.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                        help = "the path to the output file")

    filter_parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    filter_parser.add_argument("--grc", default = False, action = 'store_true',
                        help = "Use Genome Reference Consortium nomenclature rather than UCSC (default: %(default)s)")

    filter_parser.add_argument("--max_minus_log_fisher_pvalue", default = 1.0, type = float,
                       help = "remove if the - log(fisher p-value) is smaller than this value (default: %(default)s)")

    filter_parser.add_argument("--min_tumor_allele_freq", default = 0.07, type = float,
                       help = "remove if the tumor allele frequency is small than this value (default: %(default)s)")

    filter_parser.add_argument("--max_control_allele_freq", default = 0.05, type = int,
                       help = "remove if the normal allele frequency is larger than this value (default: %(default)s)")

    filter_parser.add_argument("--max_control_variant_read_pair", default = 1, type = int,
                       help = "remove if the number of variant read pairs in the matched control sample exceeds this value (default: %(default)s)")

    filter_parser.add_argument("--control_depth_thres", default = 10, type = float,
                       help = "remove if the normal read depth is smaller than this value (default: %(default)s)")

    filter_parser.add_argument("--min_overhang_size", default = 50, type = int,
                       help = "remove if either of overhang sizes for two breakpoints is below this value (default: %(default)s)")

    filter_parser.add_argument("--inversion_size_thres", default = 500, type = int,
                       help = "remove if the size of inversion is smaller than this value (default: %(default)s)")

    filter_parser.add_argument("--max_variant_size", default = None, type = int,
                       help = "remove if the size of variant is larger than this value (default: %(default)s)")

    # filter_parser.add_argument("--within_exon", default = False, action = "store_true",
    #                     help = "keep only variants within exon (default: %(default)s)")

    filter_parser.add_argument("--pooled_control_file", default = None, type = str,
                        help = "the path to control data created by merge_control (default: %(default)s)")

    filter_parser.add_argument("--pooled_control_num_thres", default = 3, type = int,
                        help = "remove if more than specified number of control sample are found (default: %(default)s)")

    filter_parser.add_argument("--simple_repeat_file", default = None, type = str,
                        help = "remove variants with overlapping simple repeat annotation. \
                            Set the annotation file downloaded from UCSC (simpleRepeat.txt.gz) (default: %(default)s)")

    filter_parser.add_argument("--remove_rna_junction", default = False, action = "store_true",
                        help = "remove putative rna splicing junction contamination (default: %(default)s)")

    filter_parser.set_defaults(func = filter_main)


    ##########
    # annotation
    annotation_parser = subparsers.add_parser("annotation", help = "annotate structural variations")

    annotation_parser.add_argument("sv_result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                               help = "the path to genomon SV result")

    annotation_parser.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                               help = "the path to the output file")

    annotation_parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    annotation_parser.add_argument("--grc", default = False, action = 'store_true',
                        help = "Use Genome Reference Consortium nomenclature rather than UCSC (default: %(default)s)")

    annotation_parser.add_argument("--re_gene_annotation", default = False, action = "store_true",
                                   help = "gene annotaiton again (default: %(default)s)")

    annotation_parser.add_argument("--closest_exon", default = False, action = "store_true",
                                   help = "add the closest exon and distance to them (default: %(default)s)")

    annotation_parser.add_argument("--closest_coding", default = False, action = "store_true",
                                   help = "add the closest coding exon and distance to them (default: %(default)s)")

    annotation_parser.add_argument("--coding_info", default = False, action = "store_true",
                                   help = "get coding information (default: %(default)s)")

    annotation_parser.add_argument("--fusion_list", default = None, type = str,
                                   help = "add known fusion information (default: %(default)s)")

    annotation_parser.set_defaults(func = annotation_main)


    ##########
    # mutation
    mutation_parser = subparsers.add_parser("mutation",
                                          help = "Add somatic SNVs and short indels to the GenomonSV results")

    mutation_parser.add_argument("sv_result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                               help = "the path to genomon SV result")

    mutation_parser.add_argument("mutation_result_file", metavar = "genomon_mutation.result.txt", default = None, type = str,
                                 help = "the path to the genomon mutation result file")

    mutation_parser.add_argument("output_file", metavar = "output.txt", default = None, type = str,
                               help = "the path to the output file")

    mutation_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                               help = "the path to the reference genome sequence")

    mutation_parser.set_defaults(func = mutation_main)

    ##########
    # concentrate
    concentrate_parser = subparsers.add_parser("concentrate",
                                               help = "List up concentrated structural variations")

    concentrate_parser.add_argument("result_list", metavar = "result_list.txt", default = None, type = str,
                                    help = "1st column: sample IDs, 2nd column: tumor type, 3rd column: genomon SV result path")

    concentrate_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                                    help = "the path to the output file")

    concentrate_parser.add_argument("--set_count", metavar = "set_count", default = 2, type = int,
                                    help = "extract #variants is equal or more than this value within specified margin (default: %(default)s)")
     
    concentrate_parser.add_argument("--set_margin", metavar = "set_margin", default = 500, type = int,
                                    help = "extract #variants is equal or more than the specified threshould within this value (default: %(default)s)")

    concentrate_parser.set_defaults(func = concentrate_main)


    ##########

    # merge control
    merge_control = subparsers.add_parser("merge_control",
                                          help = "Merge, compress and index the lists of GenomonSV results")

    merge_control.add_argument("result_list", metavar = "result_list.txt", default = None, type = str,
                               help = "1st column: sample IDs, 2nd column: tumor type, 3rd column: genomon SV result path")

    merge_control.add_argument("output_file", metavar = "merge_control.bedpe.gz", default = None, type = str,
                               help = "the path of the output file")

    merge_control.set_defaults(func = merge_control_main)

    ##########
    # realign
    realign_parser = subparsers.add_parser("realign", help = "Realign short reads around the structural variation candidates for mainly validation purpose")

    realign_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                                help = "the path to genomon SV result")
        
    realign_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                                help = "the path to the output file")

    realign_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                                help = "the path to the reference genomoe sequence")

    realign_parser.add_argument("--tumor_bam", metavar = "tumor.bam", default = None, type = str, required=True,
                                help = "the path to the tumor bam file")

    realign_parser.add_argument("--control_bam", metavar = "control.bam", default = None, type = str,
                                help = "the path to the control bam file (optional)")  
           
    realign_parser.set_defaults(func = realign_main)

    ###########
    # primer 
    primer_parser = subparsers.add_parser("primer", help = "Generate primer sequence for mainly PCR validation")

    primer_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                                help = "the path to genomon SV result")

    primer_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                                help = "the path to the output file")

    primer_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                                help = "the path to the reference genomoe sequence")

    # contig_parser.add_argument("--length", default = 250, type = int,
    #                            help = "size of each two sequence length from breakpoints")

    primer_parser.set_defaults(func = primer_main)


    ##########
    # convert to vcf format
    format_parser = subparsers.add_parser("format", help = "Convert to vcf format for short deletions and tandem duplications")

    format_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                            help = "the path to genomon SV result")

    format_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                            help = "the path to the output file")

    format_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                               help = "the path to the reference genomoe sequence")

    format_parser.add_argument('--format', choices=['vcf'], default = 'vcf',
                        help = "the format of sv file. currently we only support vcfformat")

    format_parser.add_argument("--max_size_thres", default = 100, type = int,
                            help = "remove if the size of variant is larger than this value (default: %(default)s)")

    format_parser.set_defaults(func = format_main)

    ##########
    # get homology match size

    """
    homology_parser = subparsers.add_parser("homology",
                                            help = "get micro-homology size for each SV candidate")

    homology_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                                 help = "the path to genomon SV result")

    homology_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                                 help = "the path to the output file")

    homology_parser.add_argument("--reference", metavar = "reference.fa", default = "", type = str,
                               help = "the path to the reference genome sequence (default: %(default)s)")

    homology_parser.set_defaults(func = homology_main)
    """
    ##########

    # get nonB_DB motif distance 

    nonB_DB_parser = subparsers.add_parser("nonB_DB",
                                            help = "Get distance to nonB_DB annotated region for each SV candidate")

    nonB_DB_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                                 help = "the path to genomon SV result")

    nonB_DB_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                                 help = "the path to the output file")

    nonB_DB_parser.add_argument("--nonB_DB", metavar = "nonB_DB.bed.gz", default = None, type = str, required=True,
                                help = "the path to the nonB structure database")

    nonB_DB_parser.set_defaults(func = nonB_DB_main)


    ##########
    # check the RSS motif

    RSS_parser = subparsers.add_parser("RSS",
                                       help = "Check recombination signal sequence motif near breakpoints")

    RSS_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                             help = "the path to genomon SV result")

    RSS_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                             help = "the path to the output file")

    RSS_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                            help = "the path to the reference genomoe sequence")

    RSS_parser.add_argument("--check_size", default = 50, type = float,
                            help = "check specified number of nucleotides both upstream and downstream from the breakpoints (default: %(default)s)")

    RSS_parser.set_defaults(func = RSS_main)


    ##########
    # check the AID motif

    AID_parser = subparsers.add_parser("AID",
                                       help = "Check AID motif (CG, WGCW) near breakpoints")

    AID_parser.add_argument("result_file", metavar = "genomonSV.result.txt", default = None, type = str,
                             help = "the path to genomon SV result")

    AID_parser.add_argument("output", metavar = "output.txt", default = None, type = str,
                             help = "the path to the output file")

    AID_parser.add_argument("--reference", metavar = "reference.fa", default = None, type = str, required=True,
                            help = "the path to the reference genomoe sequence")

    AID_parser.add_argument("--check_size", default = 10, type = float,
                            help = "check specified number of nucleotides both upstream and downstream from the breakpoints (default: %(default)s)")

    AID_parser.set_defaults(func = AID_main)

    ###########

    return parser


