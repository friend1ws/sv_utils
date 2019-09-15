#! /usr/bin/env python

import unittest
import os, subprocess, tempfile, shutil, filecmp
import sv_utils

class TestNonB_DB(unittest.TestCase):

    def setUp(self):
        # prepare reference genome
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        pwd_dir = os.getcwd()

        if not os.path.exists(cur_dir + "/resource"):
            os.mkdir(cur_dir + "/resource")

        if not os.path.exists(cur_dir + "/resource/nonB_DB/nonB_DB.bed.gz"):
            os.chdir(cur_dir + "/resource")
            subprocess.check_call(["git", "clone", "https://github.com/friend1ws/nonB_DB.git"])
            os.chdir("nonB_DB")
            subprocess.check_call(["bash", "prepNonBDB.sh"])
            os.chdir(pwd_dir)

        self.parser = sv_utils.parser.create_parser()


    def tearDown(self):
        # pass
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        shutil.rmtree(cur_dir + "/resource/nonB_DB/")
 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-253J-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-253J-DNA-08.genomonSV.result.nonB_DB.txt"
        nonB_DB = cur_dir + "/resource/nonB_DB/nonB_DB.bed.gz"
        answer_file = cur_dir + "/data/nonB_DB/CCLE-253J-DNA-08.genomonSV.result.nonB_DB.txt"

        # print(output_file) 
        args = self.parser.parse_args(["nonB_DB", input_file, output_file, "--nonB_DB", nonB_DB])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

