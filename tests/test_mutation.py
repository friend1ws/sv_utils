#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import sv_utils
from .check_download import *

class TestMutation(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")
       
        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        mut_file = cur_dir + "/data/CCLE-253J-DNA-08.genomon_mutation.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.mutation.txt"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/mutation/CCLE-253J-DNA-08.genomonSV.result.mutation.txt"
 
        # print(' '.join(["sv_utils", "mutation", input_file, mut_file, output_file, "--reference", ref_genome]))
        args = self.parser.parse_args(["mutation", input_file, mut_file, output_file, "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

