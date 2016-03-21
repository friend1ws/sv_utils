# sv_utils

utility scripts for processing and summarizing Genomon-SV results

## Dependency

### Python
Python (>= 2.7), `genomonsv`, `pysam`, [`primer3-py`](http://libnano.github.io/primer3-py/index.html) (for primer) packages

### Software
tabix, bgzip, blat (for realign)

## Install

```
git clone https://github.com/friend1ws/sv_utils.git
cd sv_utils
python setup.py build
python setup.py install
```

## Commands
For detailed description on each option, please consult the help for each command

### count
count the number of each variatns (deletion, tandem_duplication, inversion and translocation) for each sample
```
sv_utils count [-h] [--inseq] result_list.txt output.txt
```
### filter
process and filter GenomonSV results
```
sv_utils filter [-h] [--fisher_thres FISHER_THRES]
                       [--tumor_freq_thres TUMOR_FREQ_THRES]
                       [--normal_freq_thres NORMAL_FREQ_THRES]
                       [--normal_depth_thres NORMAL_DEPTH_THRES]
                       [--inversion_size_thres INVERSION_SIZE_THRES]
                       [--max_size_thres MAX_SIZE_THRES] [--within_exon]
                       [--control CONTROL]
                       [--control_num_thres CONTROL_NUM_THRES]
                       [--remove_simple_repeat] [--closest_exon]
                       [--closest_coding]
                       [--mutation_result genomon_mutation.result.txt]
                       [--reference reference.fa] [--re_annotation]
                       [--coding_info]
                       genomonSV.result.txt output.txt annotation_dir
```
### gene_summary
summarize the number of variants for each gene for each cancer type
```
sv_utils gene_summary [-h] [--inframe_info]
                             result_list.txt output.txt annotation_dir
                             cancer_gene_lis
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
generate primer sequeces for mainly PRC validation
```
sv_utils primer [-h] --reference reference.fa
                       genomonSV.result.txt output.txt
```


