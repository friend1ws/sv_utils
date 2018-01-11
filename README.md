# sv_utils

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/friend1ws/sv_utils.svg?branch=devel)](https://travis-ci.org/friend1ws/sv_utils)

utility scripts for processing and summarizing Genomon-SV results

## Dependency

### Python
Python (>= 2.7), [`genomonsv`](https://github.com/Genomon-Project/GenomonSV) (for merge_control, realign and primer), `pysam`, [`primer3-py`](http://libnano.github.io/primer3-py/index.html) (for primer) packages

### Software
tabix, bgzip, blat (for realign)

## Install

```
git clone https://github.com/friend1ws/sv_utils.git
cd sv_utils
pip install .
```

## Commands
For detailed description on each option, please consult the help for each command

### count
count the number of each variatns (deletion, tandem_duplication, inversion and translocation) for each sample
```
sv_utils count [-h] [--inseq] result_list.txt output.txt
```

### gene_summary
summarize the number of variants for each gene for each cancer type
```
sv_utils gene_summary [-h] [--inframe_info]
                             result_list.txt output.txt annotation_dir
                             cancer_gene_lis
```

### filter
process and filter GenomonSV results
```
sv_utils filter [-h] [--genome_id {hg19,hg38,mm10}] [--grc]
                       [--max_minus_log_fisher_pvalue MAX_MINUS_LOG_FISHER_PVALUE]
                       [--min_tumor_allele_freq MIN_TUMOR_ALLELE_FREQ]
                       [--max_control_allele_freq MAX_CONTROL_ALLELE_FREQ]
                       [--max_control_variant_read_pair MAX_CONTROL_VARIANT_READ_PAIR]
                       [--control_depth_thres CONTROL_DEPTH_THRES]
                       [--min_overhang_size MIN_OVERHANG_SIZE]
                       [--inversion_size_thres INVERSION_SIZE_THRES]
                       [--max_variant_size MAX_VARIANT_SIZE]
                       [--pooled_control_file POOLED_CONTROL_FILE]
                       [--pooled_control_num_thres POOLED_CONTROL_NUM_THRES]
                       [--simple_repeat_file SIMPLE_REPEAT_FILE]
                       [--remove_rna_junction]
                       genomonSV.result.txt output.txt
```

### annotation
Annotate GenomonSV results
```
sv_utils annotation [-h] [--genome_id {hg19,hg38,mm10}] [--grc]
                           [--re_gene_annotation] [--closest_exon]
                           [--closest_coding] [--coding_info]
                           [--fusion_list FUSION_LIST]
                           genomonSV.result.txt output.txt
```

### concentrate
obtain a list of concentrated variants
```
sv_utils concentrate [-h] [--set_count set_count]
                            [--set_margin set_margin]
                            result_list.txt output.txt
```
### merge_control
gather, compress and index control variants

```
sv_utils merge_control [-h] result_list.txt output_prefix
```
### realign
realign short reads around the structural variation candidates in the input file and check the allele frequencies
```
sv_utils realign [-h] --reference reference.fa --tumor_bam tumor.bam
                        [--control_bam control.bam]
                        genomonSV.result.txt output.txt
```

### primer
generate primer sequeces for mainly PCR validation
```
sv_utils primer [-h] --reference reference.fa
                       genomonSV.result.txt output.txt
```

### format
convert to vcf format for short deletions and tandem duplications
```
sv_utils format [-h] --reference reference.fa [--format {vcf}]
                       [--max_size_thres MAX_SIZE_THRES]
                       genomonSV.result.txt output.txt
```

### homology
get micro-homology size for each SV candidate
```
sv_utils nonB_DB [-h] genomonSV.result.txt output.txt annotation_dir
```

### nonB_DB
get nonB_DB distance for each SV candidate
```
sv_utils nonB_DB [-h] genomonSV.result.txt output.txt annotation_dir
```

### RSS
check recombination signal sequence motif near breakpoints
```
sv_utils RSS [-h] --reference reference.fa genomonSV.result.txt output.txt
```
