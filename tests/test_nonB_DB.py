#! /usr/bin/env python

import unittest
import os, subprocess, tempfile, shutil, filecmp
import sv_utils

class TestMutation(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        if not os.path.exists(cur_dir + "/nonB_DB"):
            os.mkdir(cur_dir + "/nonB_DB")

        if not os.path.exists(cur_dir + "/nonB_DB/nonB_DB.bed.gz") or not os.path.exists(cur_dir + "/nonB_DB/nonB_DB.bed.gz.tbi"):
            os.chdir(cur_dir + "/nonB_DB")
            subprocess.check_call(["bash", "prepNonBDB.sh"])
            os.chdir(cur_dir)

        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.nonB_DB.txt"
        nonB_DB = cur_dir + "/nonB_DB/nonB_DB.bed.gz"
        answer_file = cur_dir + "/data/nonB_DB/CCLE-253J-DNA-08.genomonSV.result.nonB_DB.txt"
 
        args = self.parser.parse_args(["nonB_DB", input_file, output_file, "--nonB_DB", nonB_DB])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

