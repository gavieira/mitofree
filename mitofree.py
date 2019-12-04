#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py


import argparse


class mitofree_input():
	'''This class holds and manipulates MitoFree's input data''' 
	def __init__(self, input_file):
		for line in sra_list: #Write a "parse input" function that stores its content in a suitable data strcuture
            		if line.startswith("#"):
                		continue
            		self.accession = line.split("\t")[0].strip()
            		self.species = line.split("\t")[1].strip()
            		self.sp_acc = "{}-{}".format(self.species, self.accession)
            		self.seed = line.split("\t")[2].strip()
            		#self.new_working_dir = "%s/%s-%s" % (base_working_dir, species.upper(), accession.upper())
            		self.name_of_sra_file = "%s.%s.sra" % (self.species,self.accession)
            		self.name_of_fastq_file = "%s.%s.fastq" % (self.species,self.accession)
            		self.name_of_config_file = "%s.%s.config" % (self.species,self.accession)
            		self.name_of_seed_file = "%s.seed.fa" % (self.seed)
            		self.name_of_novop_assembly_circular = "Circularized_assembly_1_%s-%s.fasta" % (self.species, self.accession) # The "1" should be changed to regex in order to accept any digit, but os.path.isfile (used in the "merge contigs" section) does not work with regex 
            		self.name_of_novop_assembly_merged = "Option_1_%s-%s.fasta" % (self.species, self.accession) #In this case, NOVOPlasty managed to merge the contigs, and if this file contains only one contig, we are going to use it for the next steps without the use of CAP3.
            		self.name_of_novop_assembly_partial = "Contigs_1_%s-%s.fasta" % (species, accession) #Partial assemblies, unmerged


def main_function(sra_list):
    with open(sra_list) as sra_list:
        base_working_dir = os.getcwd()
        print("Base working directory is '%s'" % (base_working_dir))
        for line in sra_list: #Write a "parse input" function that stores its content in a suitable data strcuture
            if line.startswith("#"):
                continue
            accession = line.split("\t")[0].strip() #e.g. "Atta_laevigata". Readable, but confusing if more than one sample from the same species are being used
            species = line.split("\t")[1].strip() #e.g. "SRR389145". Not very readable, but can be useful when using more than one sample per species
            #species_and_accession = "{}-{}".format(species, accession) #e.g. "Atta_laevigata-SRR38914"; All the advantages of species (readability) and accession (specificity)
            seed = line.split("\t")[2].strip()
            new_working_dir = "%s/%s-%s" % (base_working_dir, species.upper(), accession.upper())
            name_of_sra_file = "%s.%s.sra" % (species,accession)
            name_of_fastq_file = "%s.%s.fastq" % (species,accession)
            name_of_config_file = "%s.%s.config" % (species,accession)
            name_of_seed_file = "%s.seed.fa" % (seed)
            name_of_novop_assembly_circular = "Circularized_assembly_1_%s-%s.fasta" % (species, accession) # The "1" should be changed to regexin order to accept any digit, but os.path.isfile (used in the "merge contigs" section) does not work with regex 
            name_of_novop_assembly_merged = "Option_1_%s-%s.fasta" % (species, accession) #In this case, NOVOPlasty managed to merge the contigs, and if this file contains only one contig, we are going to use it for the next steps without the use of CAP3.
            name_of_novop_assembly_partial = "Contigs_1_%s-%s.fasta" % (species, accession) #Partial assemblies, unmerged



def getArgs():
	parser = argparse.ArgumentParser(description="Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty and MITObim")
	parser.add_argument("-S", "--savespace", action="store_true", default=False, help="Automatically removes residual assembly files such as fastq and mitobim iterations")
	parser.add_argument("-M", "--maxmemory", type=int, metavar="", default=0, help="Limit of RAM usage for NOVOPlasty. Default: no limit")
	parser.add_argument("-K", "--kmer", type=int, metavar="", default=39, help="K-mer used in NOVOPlasty assembly. Default: 39")
	parser.add_argument("-s", "--subset", type=int, metavar="", default=50000000, help="Max number of reads used in the assembly process. Default: 50 million reads")
	#parser.add_argument("-P", "--parallel", type=int, metavar="", default=1, help="Number of parallel assemblies. Default: 1")
	parser.add_argument("-T", "--timeout", type=int, metavar="", default=24, help="Custom timeout for MITObim, in hours. Default: 24h")
	parser.add_argument("filename", type=str, metavar="FILENAME", help="Path to file with multiple accessions (one per line)")
	return parser.parse_args()

if __name__ == "__main__":
	args = getArgs()

