#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess, os
from mitoattributes import mitofree_attributes


class mitoannotation(mitofree_attributes):
    def __init__(self, dataset_line, gencode=2):
        super().__init__(dataset_line)
        self.gencode = gencode
        self.refdir = "{}/refseq81m".format(self.scriptdir)


#runmitos.py -i <fasta_file> -c <genetic_code> -o <output_dir> -r <reference_dir>
    
    def run_mitos(self):
        print("Running MITOS...")
        #os.mkdir('./MITOS/')
        print("Finished MITOS with exit status {}".format(str(mitos.returncode)))
    
    
    def generate_genbank(self):
        pass