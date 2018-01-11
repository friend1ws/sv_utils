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
Summarize the frequency of each variant type (deletion, tandem_duplication, inversion and translocation) for each sample
```
sv_utils count [-h] [--inseq] result_list.txt output.txt
```

### gene_summary
Summarize the frequency of each variant type (deletion, tandem_duplication, inversion, translocation) for each cancer gene
```
gene_summary [-h] --cancer_gene_list cancer_gene_list
                             [--inframe_info]
                             result_list.txt output.txt
```
For the **--cancer_gene_list** argument, the `CancerGeneSummary/CancerGeneSummary.proc.txt` in the [`cancer_gene_db`](https://github.com/friend1ws/cancer_gene_db) repository or its variation could be used.

### filter
Filter out GenomonSV results outside specified conditions
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
For the **--fusion_list** argument, the `fusion_db.txt` file in the [`fusion_db`](https://github.com/friend1ws/fusion_db) repository cound be used.

### mutation
Add somatic SNVs and short indels to the GenomonSV results
```
sv_utils mutation [-h] --reference reference.fa
                         genomonSV.result.txt genomon_mutation.result.txt
                         output.txt
```

### concentrate
List up concentrated structural variations
```
sv_utils concentrate [-h] [--set_count set_count]
                            [--set_margin set_margin]
                            result_list.txt output.txt
```

### merge_control
Merge, compress and index the lists of GenomonSV results
```
sv_utils merge_control [-h] result_list.txt merge_control.bedpe.gz
```

### realign
Realign short reads around the structural variation candidates for mainly validation purpose
```
sv_utils realign [-h] --reference reference.fa --tumor_bam tumor.bam
                        [--control_bam control.bam]
                        genomonSV.result.txt output.txt
```

### primer
Generate primer sequence for mainly PCR validation
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


### nonB_DB
Get distance to nonB_DB annotated region for each SV candidate
```
sv_utils nonB_DB [-h] --nonB_DB nonB_DB.bed.gz
                        genomonSV.result.txt output.txt
```
For the **--nonB_DB** argument, the `nonB_DB.bed.gz` in the [`nonB_DB`](https://github.com/friend1ws/nonB_DB) repository could be used.

### RSS
Check recombination signal sequence motif near breakpoints
```
sv_utils RSS [-h] --reference reference.fa [--check_size CHECK_SIZE]
                    genomonSV.result.txt output.txt
```

### AID
Check AID motif (CG, WGCW) near breakpoints
```
sv_utils AID [-h] --reference reference.fa [--check_size CHECK_SIZE]
                    genomonSV.result.txt output.txt
```
