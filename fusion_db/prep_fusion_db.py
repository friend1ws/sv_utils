#! /usr/bin/env/python

import xlrd, re

fusion_pair = {}
with open("Mitelman/molclingene.dat", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        # print line.rstrip('\n')
        genes = sorted(F[5].split("/"))
        if len(genes) == 2:
            fusion_pair[genes[0] + '\t' + genes[1]] = "mitelman" 

with open("CosmicFusion/CosmicFusionExport.tsv", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[11] == "": continue
        if F[0] == "1822547": continue # add hoc treatment
        gene_infos = F[11].split(':')

        if len(gene_infos) < 2: continue

        genes_tmp = []
        for i in range(len(gene_infos)):
            match = re.search(r'\{[\_\w\.]+\}$', gene_infos[i])
            if match is None: continue

            gene_tmp_info = re.sub(r'\{[\_\w\.]+\}$', '', gene_infos[i])
            gene_tmp_split = gene_tmp_info.split('_')
            genes_tmp.append(gene_tmp_split[len(gene_tmp_split) - 1])

        genes_tmp = sorted(list(set(genes_tmp)))
        if len(genes_tmp) < 2: continue

        genes = sorted([genes_tmp[0], genes_tmp[1]])
        if len(genes) == 2:
            if genes[0] + '\t' + genes[1] not in fusion_pair:
                fusion_pair[genes[0] + '\t' + genes[1]] = "COSMIC"
            else:
                fusion_pair[genes[0] + '\t' + genes[1]] = fusion_pair[genes[0] + '\t' + genes[1]] + ";" + "COSMIC"

with xlrd.open_workbook("YoshiharaEtAl2015/onc2014406x2.xls") as book:
    sheet1 = book.sheet_by_index(1)
    for i in range(2, sheet1.nrows):
        genes = sorted([sheet1.cell(i, 2).value, sheet1.cell(i, 3).value])
        if genes[0] + '\t' + genes[1] not in fusion_pair:
            fusion_pair[genes[0] + '\t' + genes[1]] = "Yoshihara"
        else:
            fusion_pair[genes[0] + '\t' + genes[1]] = fusion_pair[genes[0] + '\t' + genes[1]] + ";" + "Yoshihara"


with xlrd.open_workbook("KlijnEtAl2015/nbt.3080-S7.xls") as book:
    sheet1 = book.sheet_by_index(0)
    for i in range(2, sheet1.nrows):
        genes = sorted([sheet1.cell(i, 1).value, sheet1.cell(i, 2).value])
        if genes[0] + '\t' + genes[1] not in fusion_pair:
            fusion_pair[genes[0] + '\t' + genes[1]] = "Klijn_cellline"
        else:
            fusion_pair[genes[0] + '\t' + genes[1]] = fusion_pair[genes[0] + '\t' + genes[1]] + ";" + "Klijn_cellline" 
 

with xlrd.open_workbook("KlijnEtAl2015/nbt.3080-S8.xls") as book:
    sheet1 = book.sheet_by_index(0)
    for i in range(2, sheet1.nrows):
        genes = sorted(sheet1.cell(i, 2).value.split(' '))
        if genes[0] + '\t' + genes[1] not in fusion_pair:
            fusion_pair[genes[0] + '\t' + genes[1]] = "Klijn_tcga"
        else:
            fusion_pair[genes[0] + '\t' + genes[1]] = fusion_pair[genes[0] + '\t' + genes[1]] + ";" + "Klijn_tcga"


with xlrd.open_workbook("StranskyEtAl2014/ncomms5846-s3.xlsx") as book:
    sheet1 = book.sheet_by_index(0)
    for i in range(2, sheet1.nrows):
        genes = sorted([sheet1.cell(i, 0).value, sheet1.cell(i, 1).value])
        if genes[0] + '\t' + genes[1] not in fusion_pair:
            fusion_pair[genes[0] + '\t' + genes[1]] = "Stransky"
        else:
            fusion_pair[genes[0] + '\t' + genes[1]] = fusion_pair[genes[0] + '\t' + genes[1]] + ";" + "Stransky"


for pair in sorted(fusion_pair):
    infos = ';'.join(sorted(list(set(fusion_pair[pair].split(';')))))
    print pair + '\t' + infos


