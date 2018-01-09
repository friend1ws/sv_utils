#! /usr/bin/env python

import sys, glob, os

allfiles = sorted(glob.glob("*/tsv/chr*.tsv"))

for file in allfiles:
    with open(file, 'r') as hin:
        header = hin.readline()

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "chr": continue
            print '\t'.join([F[0], F[3], F[4], F[1], F[2]] + F[5:])


