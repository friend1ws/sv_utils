#! /usr/bin/env python

import sys, gzip

sym2exist = {}
with gzip.open("../../resource/refGene.txt.gz", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        sym2exist[F[12]] = 1


hout = open("Kinase_com_manual.proc.txt", 'w')
with open("Kinase_com_manual.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        genes = [F[0]]
        alias = F[3].strip("\"").split(',')
        for gene in alias:
            genes.append(gene.strip(' '))

        sym = ""
        for gene in genes:
            if gene in sym2exist: sym = gene

        if sym == "":
            print >> sys.stderr, "No corresponding gene symbol in RefGenes: " + '\t'.join(F)
            continue


        print >> hout, sym + '\t' + ':'.join([x.strip( ' ') for x in F[2].split(":")])

 
hout.close()

