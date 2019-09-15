#! /usr/bin/env python

from __future__ import print_function
import unittest
import os, glob, tempfile, shutil, filecmp
import sv_utils

class TestMergeControl(unittest.TestCase):

    def setUp(self):
        self.parser = sv_utils.parser.create_parser()

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        all_sv_file = glob.glob(cur_dir + "/data/*.genomonSV.result.txt")
        sample2type = {}
        with open(cur_dir + "/data/sample_type.txt", 'r') as hin:
            for line in hin:
                F = line.rstrip('\n').split('\t')
                sample2type[F[0]] = F[1]
    
        with open(tmp_dir + "/CCLE.genomonSV.reulst_list.txt", 'w') as hout:
            for sv_file in sorted(all_sv_file):
                sample = os.path.basename(sv_file).replace(".genomonSV.result.txt", '')
                print("%s\t%s\t%s" % (sample, sample2type[sample], sv_file), file = hout) 

        result_lsit_file = tmp_dir + "/CCLE.genomonSV.reulst_list.txt"
        output_file = tmp_dir + "/merge_control.bedpe.gz"
        answer_file = cur_dir + "/data/merge_control/merge_control.bedpe.gz"

        args = self.parser.parse_args(["merge_control", result_lsit_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

