#! /usr/bin/env python

import os


def summary_main(args):

    hout = open(args.output, 'w')
    print >> hout, '\t'.join(["sample", "deletion", "tandem_duplication", "inversion", "translocation"])

    with open(args.result_list, 'r') as hin:

        for line in hin:
            sample, result_file = line.rstrip('\n').split('\t')
            if not os.path.exists(result_file):
                raise ValueError("file not exists: " + result_file)

            type2count = {"deletion": 0, "tandem_duplication": 0, "inversion": 0, "translocation": 0}
            with open(result_file, 'r') as hin2:
                for line in hin2:
                    F = line.rstrip('\n').split('\t')

                    if float(F[14]) < args.tumor_freq_thres: continue
                    if int(F[15]) + int(F[16]) < args.normal_depth_thres: continue
                    if float(F[17]) > args.normal_freq_thres: continue
                    if float(F[18]) < args.fisher_thres: continue 

                    type2count[F[7]] = type2count[F[7]] + 1
                       
            print >> hout, sample + '\t' + str(type2count["deletion"]) + '\t' + str(type2count["tandem_duplication"]) + '\t' + \
                                        str(type2count["inversion"]) + '\t' + str(type2count["translocation"])
        

 
    hout.close()

