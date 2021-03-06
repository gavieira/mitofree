#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess, os, shutil, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mitoattributes import mitofree_attributes

import functools
print = functools.partial(print, flush=True)


class mitoannotation(mitofree_attributes):
    def __init__(self, dataset_line, gencode=2):
        super().__init__(dataset_line)
        self.contigsize = len(SeqIO.read(self.mitobim_result, 'fasta').seq)
        self.gencode = gencode
        self.refdir = "{}/refseq63m".format(self.scriptdir)
        self.beddir = "./mitos_results/result.bed"
        self.feat_dict = {
        "cox1": {"name" : "COX1", "product" : "cytochrome c oxidase subunit I"},
        "cox2": {"name": "COX2", "product": "cytochrome c oxidase subunit II"},
        "atp8": {"name": "ATP8", "product": "ATP synthase F0 subunit 8"},
        "atp6": {"name": "ATP6", "product": "ATP synthase F0 subunit 6"},
        "cox3": {"name": "COX3", "product": "cytochrome c oxidase subunit III"},
        "nad3": {"name":"ND3", "product": "NADH dehydrogenase subunit 3"},
        "nad5": {"name": "ND5", "product": "NADH dehydrogenase subunit 5"},
        "nad4l": {"name": "ND4L", "product": "NADH dehydrogenase subunit 4L"},
        "nad4": {"name":"ND4", "product": "NADH dehydrogenase subunit 4"},
        "nad6": {"name":"ND6", "product":"NADH dehydrogenase subunit 6"},
        "cob": {"name":"CYTB", "product":"cytochrome b"},
        "nad1": {"name":"ND1", "product":"NADH dehydrogenase subunit 1"},
        "nad2": {"name":"ND2", "product":"NADH dehydrogenase subunit 2"},
        "trnL2": {"name":"trnL2", "product":"tRNA-Leu"},
        "trnK": {"name":"trnK", "product":"tRNA-Lys"},
        "trnK": {"name":"trnK", "product":"tRNA-Lys"},
        "trnD": {"name":"trnD", "product":"tRNA-Asp"},
        "trnG": {"name":"trnG", "product":"tRNA-Gly"},
        "trnA": {"name":"trnA", "product":"tRNA-Ala"},
        "trnR": {"name":"trnR", "product":"tRNA-Arg"},
        "trnN": {"name":"trnN", "product":"tRNA-Asn"},
        "trnS1": {"name":"trnS1", "product":"tRNA-Ser"},
        "trnE": {"name":"trnE", "product":"tRNA-Glu"},
        "trnF": {"name":"trnF", "product":"tRNA-Phe"},
        "trnH": {"name":"trnH", "product":"tRNA-His"},
        "trnT": {"name":"trnT", "product":"tRNA-Thr"},
        "trnP": {"name":"trnP", "product":"tRNA-Pro"},
        "trnS2": {"name":"trnS2", "product":"tRNA-Ser"},
        "trnL1": {"name":"trnL1", "product":"tRNA-Leu"},
        "trnV": {"name":"trnV", "product":"tRNA-Val"},
        "trnM": {"name":"trnM" , "product":"tRNA-Met"},
        "trnI": {"name":"trnI", "product":"tRNA-Ile"},
        "trnQ": {"name":"trnQ", "product":"tRNA-Gln"},
        "trnW": {"name":"trnW", "product":"tRNA-Trp"},
        "trnC": {"name":"trnC", "product":"tRNA-Cys"},
        "trnY": {"name":"trnY", "product":"tRNA-Tyr"},
        "rrnL": {"name":"rrnL", "product":"large subunit ribosomal RNA"},
        "rrnS": {"name":"rrnS", "product":"small subunit ribosomal RNA"}
        }
        self.cds_list = [x for x in self.feat_dict.keys() if not x.startswith("trn") and not x.startswith("rrn")]
        self.trna_list = [x for x in self.feat_dict.keys() if x.startswith("trn")]
        self.rrna_list = [x for x in self.feat_dict.keys() if x.startswith("rrn")]
        self.gbk = "{}.gbk".format(self.prefix)
        
    
    def check_mitos_results(self):
        if os.path.isfile(self.beddir):
            return True
        else:
            if os.path.isdir("mitos_results"):
                shutil.rmtree("mitos_results")
            os.mkdir("mitos_results")
            return False

    def run_mitos(self):
        if not self.check_mitos_results():
            print("Running MITOS...")
            with open("mitos.out", "w+") as output, open('mitos.err', 'w+') as error:
                mitos = subprocess.run(["runmitos.py", "-i", self.mitobim_result, "-c", str(self.gencode), "-o", "mitos_results", "-r", self.refdir, "--linear", "--ncbicode", "--best", "--alarab", "--intron", "0", "--oril", "0", "--orih", "0"], stdout=output, stderr=error)
                print("Finished MITOS with exit status {}".format(str(mitos.returncode)))
                output.seek(0)
                missing_feat = False
                for line in output:
                    if line.startswith("missing:"):
                        missing_feat = True
                        print('WARNING: {}'.format(line))
                if not missing_feat:
                    print("All features were succesfully annotated")
        else:
            print("Annotation process already finished. Skipping to generation of genbank file...")

    def generate_gbk(self):
        '''Generates a gbk file in the assembly directory'''
        with open(self.gbk, "w") as gbk:
            gbk.write(self.format_gbk_header() + self.format_features() + self.format_sequence())
        

    def format_features(self):
        formatted_feats = ''
        with open(self.beddir) as bed:
            for feature in bed:
                (feature_name, feature_type, product, anticodon, inipos, endpos, strand) = self.parse_bed(feature)
                inipos = int(inipos) + 1
                if strand.strip() == "+":
                    formatted_feats += "{}{:<16}{:<}\n".format(5*" ", feature_type, "{}..{}".format(str(inipos), endpos))
                if strand.strip() == "-":
                    formatted_feats += "{}{:<16}{:<}\n".format(5*" ", feature_type, "complement({}..{})".format(str(inipos), endpos))
                if anticodon:
                    formatted_feats += "{}{:<}\n".format(21*" ", '/note="anticodon:{}"'.format(anticodon))
                formatted_feats += "{}{:<}\n".format(21*" ", '/product="{}"'.format(product))
                formatted_feats += "{}{:<}\n".format(21*" ", '/gene="{}"'.format(feature_name))
        return formatted_feats

    def format_gbk_header(self):
        gbk_header_template = "{}/gbk_header_template.txt".format(self.scriptdir)
        with open(gbk_header_template) as gbk_header: ##OBS: The newline in the gbk_header_template.txt is important!!!
            return gbk_header.read().format(self.sra_run_number, self.contigsize, self.species)

    
    def format_sequence(self):
        formatted_seq = 'ORIGIN\n'
        mitoseq = SeqIO.read(self.mitobim_result, "fasta").seq.lower()
        subsequences = ((mitoseq[0+i:60+i], i) for i in range(0, len(mitoseq), 60))
        for subseq, index in subsequences:
            formatted_seq += "{:>10} {}\n".format(index+1, subseq)
        formatted_seq += "//"
        return formatted_seq
    
    def parse_bed(self, feature):
        line = feature.split("\t")
        (feature_name, feature_type, product, anticodon) = self.parse_feat_name(line[3])
        inipos = line[1]
        endpos = line[2]
        strand = line[5]
        return (feature_name, feature_type, product, anticodon, inipos, endpos, strand)
        
    def parse_feat_name(self, feat_name):
        anticodon = final_feat_name = feat_type = product = ''
        if feat_name.strip().startswith("trn"): #'(' index used to separate codon/name - e.g. "trnL1(tag)" (MITOS output in .bed file) to "trnL1" and "tag". The codon is then reverse translated to obtain the anticodon sequence ('cua', in this case).
            separator_index = feat_name.find("(")
            codon = feat_name[separator_index+1:-1]
            anticodon = str(Seq(codon, generic_dna).reverse_complement().transcribe())
            feat_base_name = feat_name[:separator_index]
            final_feat_name = "{}-{}".format(feat_base_name, anticodon)
            product = self.feat_dict.get(feat_base_name).get("product")
            feat_type = 'tRNA' 
        elif feat_name.startswith("rrn"):
            final_feat_name = feat_name
            product = self.feat_dict.get(feat_name).get("product")
            feat_type = 'rRNA'
        elif feat_name in self.cds_list:
            final_feat_name = self.feat_dict.get(feat_name).get("name")
            product = self.feat_dict.get(feat_name).get("product")
            feat_type = 'CDS'
        return (final_feat_name, feat_type, product, anticodon)
    
    def copy_gbk_to_new_directory(self):
        #if self.gbk is not missing any features:
        gbk_dir = "../gbk_files"
        if not os.path.exists(gbk_dir):
            os.mkdir(gbk_dir)
        shutil.copy(self.gbk, gbk_dir)
        pass        

import argparse, traceback

def getArgs():
    parser = argparse.ArgumentParser(description="Anotates mitochondrial contigs")
    parser.add_argument("filename", type=str, metavar="FILENAME", help="Path to dataset list")
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    with open(args.filename) as datasets:
        for line in datasets:
            if line.startswith('#'): pass
            else:
                try:
                    annotation = mitoannotation(line, gencode=5)
                    annotation.run_mitos()
                    annotation.generate_gbk()
                except Exception as error:
                    fullerror = traceback.format_exc()
                    print("An error has occurred for this assembly: {}\n\nFULL ERROR:\n\n{}\n\nProceeding to the next assembly...\n".format(error, fullerror))
                    pass