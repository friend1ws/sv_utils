#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import sv_utils
from check_reference_genome import *

class TestPrimer(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_reference_genome("https://storage.googleapis.com/genomon_rna_gce/db/GRCh37/GRCh37.fa", \
                                 cur_dir + "/resource/reference_genome/GRCh37.fa")
       
        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.primer.txt"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/primer/CCLE-253J-DNA-08.genomonSV.result.primer.txt"
 
        args = self.parser.parse_args(["primer", input_file, output_file, "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

