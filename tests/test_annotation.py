#! /usr/bin/env python

import unittest
import os, urllib2, tempfile, shutil, filecmp
import sv_utils

class TestAnnotation(unittest.TestCase):

    def setUp(self):
        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-K-562-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-K-562-DNA-08.genomonSV.result.annotation.txt"
        answer_file = cur_dir + "/data/annotation/CCLE-K-562-DNA-08.genomonSV.result.annotation.txt"

        args = self.parser.parse_args(["annotation", input_file, output_file, "--grc", "--re_gene_annotation", \
                                       "--closest_exon", "--closest_coding", "--coding_info", "--fusion_info"])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

