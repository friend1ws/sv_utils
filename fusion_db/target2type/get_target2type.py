#! /usr/bin/env python

import sys, xlrd

target_gene_file = sys.argv[1]
output_file = sys.argv[2]

target_gene = {}
with open(target_gene_file, 'r') as hin:
    for line in hin:
        target_gene[line.rstrip('\n')] = 1

target_gene2type = {}

# Yoshihara et al. Oncogene, 2015
with xlrd.open_workbook("../YoshiharaEtAl2015/onc2014406x2.xls") as book:
    sheet1 = book.sheet_by_index(1)
    for i in range(2, sheet1.nrows):
        genes = sorted([sheet1.cell(i, 2).value, sheet1.cell(i, 3).value])
        if genes[0] in target_gene:
            if genes[0] not in target_gene2type: target_gene2type[genes[0]] = []
            target_gene2type[genes[0]].append(sheet1.cell(i, 0).value)

        if genes[1] in target_gene:
            if genes[1] not in target_gene2type: target_gene2type[genes[1]] = []
            target_gene2type[genes[1]].append(sheet1.cell(i, 0).value)

# Klijn et al., Nature Biotechnology, 2015
with xlrd.open_workbook("../KlijnEtAl2015/nbt.3080-S8.xls") as book:
    sheet1 = book.sheet_by_index(0)
    for i in range(2, sheet1.nrows):
        genes = sorted(sheet1.cell(i, 2).value.split(' '))
        if genes[0] in target_gene:
            if genes[0] not in target_gene2type: target_gene2type[genes[0]] = []
            target_gene2type[genes[0]].append(sheet1.cell(i, 1).value)

        if genes[1] in target_gene:
            if genes[1] not in target_gene2type: target_gene2type[genes[1]] = []
            target_gene2type[genes[1]].append(sheet1.cell(i, 1).value)


# Stransky et al., Nature Communication, 2014 
stransky_type = {
"Thyroid carcinoma:  TCGA": "TYCA",
"Lung adenocarcinoma:  TCGA": "LUAD",
"Rectum adenocarcinoma:  TCGA": "READ",
"Kidney renal papillary cell carcinoma:  TCGA": "KIRP",
"Bladder Urothelial Carcinoma:  TCGA": "BLCA",
"Skin Cutaneous Melanoma:  TCGA": "SKCM",
"Prostate adenocarcinoma:  TCGA": "PRAD",
"Stomach adenocarcinoma:  TCGA": "STAD",
"Breast invasive carcinoma:  TCGA": "BRCA",
"Lung squamous cell carcinoma:  TCGA": "LUSC",
"Ovarian serous cystadenocarcinoma:  TCGA": "OV",
"Liver hepatocellular carcinoma:  TCGA": "LIHC", 
"Brain Lower Grade Glioma:  TCGA": "LGG",
"Uterine Corpus Endometrial Carcinoma:  TCGA": "UCEC",
"Glioblastoma multiforme:  TCGA": "GBM",
"Sarcoma:  TCGA": "SARC", 
"Head and Neck squamous cell carcinoma:  TCGA": "HNSC",
"Colon adenocarcinoma:  TCGA": "COAD"
}

with xlrd.open_workbook("../StranskyEtAl2014/ncomms5846-s3.xlsx") as book:
    sheet1 = book.sheet_by_index(0)
    for i in range(2, sheet1.nrows):
        genes = sorted([sheet1.cell(i, 0).value, sheet1.cell(i, 1).value])
        if genes[0] in target_gene:
            if genes[0] not in target_gene2type: target_gene2type[genes[0]] = []
            target_gene2type[genes[0]].append(stransky_type[sheet1.cell(i, 2).value])

        if genes[1] in target_gene:
            if genes[1] not in target_gene2type: target_gene2type[genes[1]] = []
            target_gene2type[genes[1]].append(stransky_type[sheet1.cell(i, 2).value])

hout = open(output_file, 'w')
for gene in sorted(target_gene):
    types = "---"
    if gene in target_gene2type:
        types = ';'.join(sorted(list(set(target_gene2type[gene]))))
    print >> hout, gene + '\t' + types
 
hout.close()
