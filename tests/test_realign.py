#! /usr/bin/env python

import unittest
import os, tempfile, shutil, filecmp
import sv_utils
from check_download import *

class TestRSS(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        check_download("https://storage.googleapis.com/friend1ws_package_data/common/GRCh37.fa", \
                       cur_dir + "/resource/reference_genome/GRCh37.fa")

        check_download("https://storage.googleapis.com/friend1ws_package_data/sv_utils/CCLE-253J-DNA-08.markdup.bam", \
                       cur_dir + "/resource/realign_bam/CCLE-253J-DNA-08.markdup.bam")

        check_download("https://storage.googleapis.com/friend1ws_package_data/sv_utils/CCLE-253J-DNA-08.markdup.bam.bai", \
                       cur_dir + "/resource/realign_bam/CCLE-253J-DNA-08.markdup.bam.bai")
 
        self.parser = sv_utils.parser.create_parser()


    def tearDown(self):
       cur_dir = os.path.dirname(os.path.abspath(__file__))
       shutil.rmtree(cur_dir + "/resource/realign_bam")
 

    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.realign.txt"
        tumor_bam_file = cur_dir + "/resource/realign_bam/CCLE-253J-DNA-08.markdup.bam"
        ref_genome = cur_dir + "/resource/reference_genome/GRCh37.fa"
        answer_file = cur_dir + "/data/realign/CCLE-253J-DNA-08.genomonSV.result.realign.txt"

        args = self.parser.parse_args(["realign", "--reference", ref_genome, "--tumor_bam", tumor_bam_file,
                                       input_file, output_file])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

