#! /usr/bin/env python

from __future__ import print_function
import unittest
import os, subprocess, tempfile, shutil, filecmp
import sv_utils

class TestAnnotation(unittest.TestCase):

    def setUp(self):

        # download cancer_gene_db
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        pwd_dir = os.getcwd()

        if not os.path.exists(cur_dir + "/resource"):
            os.makedirs(cur_dir + "/resource")

        if not os.path.exists(cur_dir + "/resource/fusion_db/fusion_db.txt"):
            os.chdir(cur_dir + "/resource")
            subprocess.check_call(["git", "clone", "https://github.com/friend1ws/fusion_db.git"])
            os.chdir("fusion_db")
            print(os.getcwd())
            print(' '.join(["bash", "prep_fusion_db.sh"]))
            subprocess.check_call(["bash", "prep_fusion_db.sh"])
            os.chdir(pwd_dir)

        self.parser = sv_utils.parser.create_parser()

 
    def test1(self):

        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()

        input_file = cur_dir + "/data/CCLE-K-562-DNA-08.genomonSV.result.txt"
        output_file = tmp_dir + "/CCLE-K-562-DNA-08.genomonSV.result.annotation.txt"
        fusion_list_file = cur_dir + "/resource/fusion_db/fusion_db.txt"
        answer_file = cur_dir + "/data/annotation/CCLE-K-562-DNA-08.genomonSV.result.annotation.txt"

        # print(output_file)
        args = self.parser.parse_args(["annotation", input_file, output_file, "--grc", "--re_gene_annotation", \
                                       "--closest_exon", "--closest_coding", "--coding_info", "--fusion_list", fusion_list_file])
        args.func(args)

        # print(output_file)
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        
        shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    unittest.main()

