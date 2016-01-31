#! /bin/env python

import sys, gzip

refgene_file = sys.argv[1]


with gzip.open(refgene_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        
        starts = F[9].split(',')
        ends = F[10].split(',')
    
        if F[1].startswith("NM_"):

            for i in range(len(starts) - 1):
                if i >= 1:
                    key = F[2] + '\t' + ends[i - 1] + '\t' + starts[i] + '\t' + F[12] + ';' + F[1]
                    print key + '\t' + "intron" + '\t' + F[3]

                if min(int(ends[i]), int(F[6])) - int(starts[i]) > 0:
                    key = F[2] + '\t' + starts[i] + '\t' + str(min(int(ends[i]), int(F[6]))) + '\t' + F[12] + ";" + F[1]
                    if F[3] == '+':
                        print key + '\t' + "5UTR" + '\t' + '+'
                    else:
                        print key + '\t' + "3UTR" + '\t' + '-'

                if min(int(ends[i]), int(F[7])) - max(int(starts[i]), int(F[6])) > 0:
                    key = F[2] + '\t' + str(max(int(starts[i]), int(F[6]))) + '\t' + str(min(int(ends[i]), int(F[7]))) + '\t' + F[12] + ";" + F[1]
                    print key + '\t' + "coding" + '\t' + F[3]

                if int(ends[i]) - max(int(F[7]), int(starts[i])) > 0:
                    key = F[2] + '\t' + str(max(int(F[7]), int(starts[i]))) + '\t' + ends[i] + '\t' + F[12] + ';' + F[1]
                    if F[3] == '+':
                        print key + '\t' + "3UTR" + '\t' + '+'
                    else:
                        print key + '\t' + "5UTR" + '\t' + '-'


        if F[1].startswith("NR_"):
        
            for i in range(len(starts) - 1):
                if i >= 1:
                    key = F[2] + '\t' + ends[i - 1] + '\t' + starts[i] + '\t' + F[12] + ';' + F[1]
                    print key + '\t' + "intron" + '\t' + F[3]

                key = F[2] + '\t' + starts[i] + '\t' + ends[i] + '\t' + F[12] + ';' + F[1]
                print key + '\t' + "noncoding" + '\t' + F[3]

        

