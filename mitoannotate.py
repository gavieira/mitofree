#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess, os, shutil
from mitoattributes import mitofree_attributes


class mitoannotation(mitofree_attributes):
    def __init__(self, dataset_line, gencode=2):
        super().__init__(dataset_line)
        self.gencode = gencode
        self.refdir = "{}/refseq81m".format(self.scriptdir)


#runmitos.py -i <fasta_file> -c <genetic_code> -o <output_dir> -r <reference_dir>
    
    def check_mitos_results(self):
        annotation_path = "./mitos_results/result.gff"
        if os.path.isfile(annotation_path):
            return True
        else:
            if os.path.isdir("mitos_results"):
                shutil.rmtree("mitos_results")
            os.mkdir("mitos_results")
            return False

    
    def run_mitos(self):
        if not self.check_mitos_results():
            print("Running MITOS...")
            mitos = subprocess.run(["runmitos.py", "-i", self.mitobim_result, "-c", str(self.gencode), "-o", "mitos_results", "-r", self.refdir])
            print("Finished MITOS with exit status {}".format(str(mitos.returncode)))
        else:
            print("Annotation process already finished. Skipping to generation of genbank file...")

    
    
    def generate_genbank(self):
        pass