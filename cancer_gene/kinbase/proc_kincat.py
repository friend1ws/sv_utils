#! /usr/bin/env python

"""
import gzip

ref2sym = {}
with gzip.open("../../resource/refGene.txt.gz", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        ref2sym[F[1]] = F[12]
"""
with open("Kincat_Hsap.08.02.txt", 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "Name": continue
        if len(F) < 23: continue
        if F[23] != "":
            print F[23] + '\t' + F[1]
 
