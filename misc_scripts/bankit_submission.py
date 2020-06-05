#!/usr/bin/env python3

from Bio import SeqIO
import argparse, os

import warnings
from Bio import BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonParserWarning) #Ignore seqIO parser warnings regarding incomplete LOCUS line

class bankit_submission():
    def __init__(self, gbk_dir, gencode, outdir):
        self.gbk_dir = gbk_dir
        self.gencode = gencode
        self.gbk_filelist = self.get_gbk_filelist()
        self.outdir = outdir
        self.sequin = os.path.join(self.outdir, "submission.sequin")
        self.fasta = os.path.join(self.outdir, "submission.fasta")
        self.sra_table = os.path.join(self.outdir, "submission_sra_table.txt")
        self.scriptdir = os.path.dirname(os.path.realpath(__file__))
        
    def get_gbk_filelist(self):
        return [os.path.join(self.gbk_dir, file) for file in os.listdir(self.gbk_dir)
                      if file.endswith('.gbk') or file.endswith('.gb')]
    
    def create_outdir(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
    
    def remove_previous_submission(self):
        for file in [self.sequin, self.fasta, self.sra_table]:
            if os.path.exists(file):
                os.remove(file)

    def format_sequin(self, gbk_seqio, seqid, organism):
        formatted_sequin = ''
        with open("{}/sequin_minimum_template.txt".format(self.scriptdir)) as sequin_template:
            content = sequin_template.readlines()
            formatted_sequin += content[0].format(seqid, organism)
            for feature in gbk_seqio.features:
                if feature.type in ["gene", "source"]:
                    continue
                location_start = feature.location.start.position + 1 #SeqIO.read always deduces one nt from the feature's start position: https://github.com/biopython/biopython/issues/897
                gene_template = content[1] + content[2]
                gene_feature = gene_template.format(location_start, feature.location.end,
                                                    "gene", 
                                                    feature.qualifiers.get("gene")[0])
                standard_template = content[1] + content[2] + content[3]
                standard_feature = standard_template.format(location_start, feature.location.end,
                                                    feature.type, 
                                                    feature.qualifiers.get("gene")[0], 
                                                    feature.qualifiers.get("product")[0])
                cds_fields = content[4] + content[5].format(self.gencode)
                if feature.type in ["tRNA", "rRNA"]:
                    formatted_sequin += standard_feature
                if feature.type == "CDS":
                    formatted_sequin += gene_feature + standard_feature + cds_fields
                if feature.qualifiers.get("note"):
                    formatted_sequin += content[6].format(feature.qualifiers.get("note")[0])
        return formatted_sequin
  
    def generate_sequin(self, formatted_sequin):
        with open(self.sequin, 'a') as sequin_out:
            sequin_out.write(formatted_sequin)

    def generate_fasta(self, gbk_seqio, seqid, organism):
        gbk_seqio.id = "{} [organism={}]".format(seqid, organism)
        with open(self.fasta, 'a') as fasta_out:
            fasta_out.write(gbk_seqio.format('fasta'))
    
    def create_SRA_table_header(self):
        with open(self.sra_table, 'w') as sra_table:
            sra_table.write("Sequence_ID\tAccession\n")
                            
    def generate_SRA_table(self, seqid, sra_accession):
        '''Only useful in TPA (Third Party Annotation) submissions: Generates file with sequence id and its corresponding SRA accession
        Note that this accession needs to be in the 'DBLINK' field of the genbank files
        If no sra_accession is found on this field, it is replaced by the string "NO_SRA_ACCESSION"'''
        with open(self.sra_table, 'a') as sra_table:
            sra_table.write("{}\t{}\n".format(seqid, sra_accession))
                        
    def prepare_submission(self):
        print("Preparing bankit submission...")
        self.remove_previous_submission()
        self.create_outdir()
        self.create_SRA_table_header()
        for counter, gbk in enumerate(self.gbk_filelist, 1):
            seqid = "Seq{}".format(str(counter))
            gbk_seqio = SeqIO.read(gbk, "genbank")
            organism = gbk_seqio.annotations.get("organism").replace("_", " ").replace("-", " ")
            formatted_sequin = self.format_sequin(gbk_seqio, seqid, organism)
            self.generate_sequin(formatted_sequin)
            self.generate_fasta(gbk_seqio, seqid, organism)
            sra_accession = "NO_SRA_ACCESSION" #We assume at first that there is no sra accession...
            if gbk_seqio.dbxrefs:
                sra_accession = gbk_seqio.dbxrefs[0] #Except when there actually is one!
            self.generate_SRA_table(seqid, sra_accession)

def getArgs():
    parser = argparse.ArgumentParser(description="Generates bankit submission files for all gbk files (.gb or .gbk) in a given directory")
    parser.add_argument("-o", "--outdir", type=str, metavar="", default="./bankit_submission", help="Output directory. Default: ./bankit_submission")
    parser.add_argument("-g", "--gencode", type=int, metavar="", default=2, help="Genetic code table. Default: 2 (Vertebrate Mitochondrial)")
    parser.add_argument("gbk_dir", type=str, metavar="gbk_dir", help="Path to directory containing genbank file(s)")
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    sub = bankit_submission(args.gbk_dir, args.gencode, args.outdir)
    sub.prepare_submission()
