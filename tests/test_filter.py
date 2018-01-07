#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import sv_utils

class TestFilter(unittest.TestCase):

    def setUp(self):
        self.parser = sv_utils.parser.create_parser()

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.filt1.txt"
        answer_file = cur_dir + "/data/filter/CCLE-253J-DNA-08.genomonSV.result.filt1.txt"
 
        args = self.parser.parse_args(["filter", input_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test2(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.filt2.txt"
        answer_file = cur_dir + "/data/filter/CCLE-253J-DNA-08.genomonSV.result.filt2.txt"
        simple_repeat_file = cur_dir + "/data/filter/simpleRepeat.light.txt.gz"

        args = self.parser.parse_args(["filter", input_file, output_file, \
                                       "--grc", "--simple_repeat_file", simple_repeat_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


    def test3(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.filt3.txt"
        answer_file = cur_dir + "/data/filter/CCLE-253J-DNA-08.genomonSV.result.filt3.txt"
        simple_repeat_file = cur_dir + "/data/filter/simpleRepeat.light.txt.gz"

        args = self.parser.parse_args(["filter", input_file, output_file, \
                                       "--grc", "--remove_rna_junction", "--simple_repeat_file", simple_repeat_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

