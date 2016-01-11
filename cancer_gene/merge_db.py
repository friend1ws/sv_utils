#! /usr/bin/env python

gene2lawrence = {}
with open("LawrenceEtAl_sig_gene.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2lawrence[F[0]] = F[1]

gene2cgc = {}
with open("cancer_gene_census.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        gene2cgc[F[0]] = F[1].strip('"')


genes = list(set(gene2lawrence.keys() + gene2cgc.keys()))


hout = open("cancer_gene.txt", 'w')
for gene in sorted(genes):
    lawrence = gene2lawrence[gene] if gene in gene2lawrence else "---"
    cgc = gene2cgc[gene] if gene in gene2cgc else "---"
    print >> hout, gene + '\t' + lawrence + '\t' + cgc

hout.close()

