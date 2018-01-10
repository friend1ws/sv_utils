#! /usr/bin/env python

import unittest
import os, subprocess, glob, tempfile, shutil, filecmp
import sv_utils

class TestGeneSummary(unittest.TestCase):

    def setUp(self):

        # download cancer_gene_db
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        pwd_dir = os.getcwd()
 
        if not os.path.exists(cur_dir + "/resource/cancer_gene_db/CancerGeneSummary/CancerGeneSummary.proc.txt"):
            os.chdir(cur_dir + "/resource")
            subprocess.check_call(["git", "clone", "https://github.com/friend1ws/cancer_gene_db.git"])
            os.chdir("cancer_gene_db")
            subprocess.check_call(["bash", "prep_all.sh"])
            os.chdir(pwd_dir)

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
                print >> hout, "%s\t%s\t%s" % (sample, sample2type[sample], sv_file) 

        result_lsit_file = tmp_dir + "/CCLE.genomonSV.reulst_list.txt"
        cancer_gene_list = cur_dir + "/resource/cancer_gene_db/CancerGeneSummary/CancerGeneSummary.proc.txt"
        output_file = tmp_dir + "/gene_summary_sv_list.txt"
        answer_file = cur_dir + "/data/gene_summary/gene_summary_sv_list.txt"


        args = self.parser.parse_args(["gene_summary", result_lsit_file, output_file, "--cancer_gene_list", cancer_gene_list])
        args.func(args)

        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))

        shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    unittest.main()

