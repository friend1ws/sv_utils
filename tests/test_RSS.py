#! /usr/bin/env python

import unittest
import os, urllib2, tempfile, shutil, filecmp
import sv_utils

class TestRSS(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(cur_dir + "/reference_genome"):
            os.mkdir(cur_dir + "/reference_genome")

        if not os.path.exists(cur_dir + "/reference_genome/GRCh37.fa"):
            ufile = urllib2.urlopen('https://storage.googleapis.com/genomon_rna_gce/db/GRCh37/GRCh37.fa')
            with open(cur_dir + "/reference_genome/GRCh37.fa",'w') as hout:
                for x in ufile:
                    hout.write(x)
       
        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.RSS.txt"
        ref_genome = cur_dir + "/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/RSS/CCLE-253J-DNA-08.genomonSV.result.RSS.txt"
 
        args = self.parser.parse_args(["RSS", input_file, output_file, "--reference", ref_genome])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

