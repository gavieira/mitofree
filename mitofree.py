#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"


#accession = SRR5437752 (length 10)
#ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR543/002/SRR5437752
#Pattern = ftp://ftp.sra.ebi.ac.uk/vol1/accession[0:3].lower()/accession[0:6]/00accession[-1]/accession

#accession = ERR020102 (length 9)
#ftp://ftp.sra.ebi.ac.uk/vol1/err/ERR969/ERR969522
#Pattern = ftp://ftp.sra.ebi.ac.uk/vol1/accession[0:3].lower()/accession[0:6]//accession

import argparse
parser = argparse.ArgumentParser(description="Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty and MITObim")
parser.add_argument("-M", "--maxmemory", type=int, metavar="", default=0, help="Limit of RAM usage for NOVOPlasty. Default: no limit")
parser.add_argument("-K", "--kmer", type=int, metavar="", default=39, help="K-mer used in NOVOPlasty assembly. Default: 39")
parser.add_argument("filename", type=str, metavar="FILENAME", help="Path to file with multiple accessions (one per line)")
args = parser.parse_args()

import wget
import os
import re
from Bio import Entrez

def main_function(sra_list):
    with open(sra_list) as sra_list:
        base_working_dir = os.getcwd()
        print("Base working directory is '%s'" % (base_working_dir))
        for line in sra_list:    
            accession = line.split("\t")[0].strip()
            species = line.split("\t")[1].strip()
            seed = line.split("\t")[2].strip()
            new_working_dir = "%s/%s-%s" % (base_working_dir, species.upper(), accession.upper())
            name_of_sra_file = "%s.%s.sra" % (species,accession)
            name_of_fastq_file = "%s.%s.fastq" % (species,accession)
            name_of_config_file = "%s.%s.config" % (species,accession)
            name_of_seed_file = "%s.seed.fa" % (seed)
            if create_folders(new_working_dir):
                if download_sra_files(accession, name_of_sra_file, new_working_dir):
                    max_read_length = highest_read_length(name_of_sra_file, name_of_fastq_file)
                    generate_fastq(name_of_sra_file, max_read_length)
                    download_seed(name_of_seed_file, seed)
                    run_NOVOPlasty(accession, species, name_of_fastq_file, name_of_config_file, name_of_seed_file,max_read_length)
        return("All done!")


def generate_ftp_link(accession):
    url = ""
    if len(accession) == 10:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/%s/%s/00%s/%s" % (accession[:3].lower(), accession[:6], accession[-1], accession)
    if len(accession) == 9:
        url = "ftp://ftp.sra.ebi.ac.uk/vol1/%s/%s/%s" % (accession[:3].lower(), accession[:6], accession)
    return url

def create_folders(new_working_dir): ##Creates folder for each dataset and changes the working directory
    print("New working directory is '%s'\n" % (new_working_dir))
    try:
        output_mkdir = os.system("mkdir -p %s" % (new_working_dir))
        os.chdir(new_working_dir)
        return True
    except re.search(".*Permission denied$", output_mkdir):
        print("Could not create folder %s. Permission denied." % (new_working_dir.split("/")[-1]))
        return False
    except:
        print("Could not create folder %s" % (new_working_dir.split("/")[-1]))
        return False

def download_sra_files(accession, name_of_sra_file, new_working_dir):
    if os.path.isfile(name_of_sra_file): ##Checks if the sra file has already been downloaded.
        print("The file %s has already been downloaded.\n" %(name_of_sra_file))
        return False
    else:
        try:
            print("Downloading %s:" % (accession))
            wget.download(generate_ftp_link(accession), out= name_of_sra_file) ##Downloads sra files
            print("\n")
            return True
        except:
            os.system("rm -r %s" % (new_working_dir)) ##There is no need for the folder if no dataset has been downloaded 
            print("The %s dataset could not be downloaded. Is the run number correct?\n" % (accession))
            return False


def download_seed(name_of_seed_file, seed):
    try:
        with open(name_of_seed_file, "w+") as fasta:
            handle = Entrez.efetch(db='nucleotide', id=seed, rettype='fasta', retmode='text')
            fasta.write(handle.read())
        print("Seed file %s downloaded succesfully" % (name_of_seed_file))
    except:
        print("Seed file %s could not be downloaded" % (name_of_seed_file))

def generate_fastq(name_of_sra_file, max_read_length):
    print("Converting %s to fastq..." % (name_of_sra_file))
    try:
        os.system("fastq-dump -M %d --split-spot --defline-seq '@$ac-$sn/$ri' --defline-qual '+' -O ./ %s" % (max_read_length-1,name_of_sra_file))
        print("Dataset has been converted to fastq succesfully!\n")
    #fastq_name = re.sub("sra$", "fastq", sra_file)
    #with open(fastq_name, "a") as fastq:
    except:
        print("Dataset could not be converted to fastq\n")
    #use sratoolkit's fastq-dump to generate a subset(or not) of reads with header
    #use regex to extract list of read sizes "length=([0-9]+)$"
    #identify max read length and use it as parameter for fastq-dump


def highest_read_length(name_of_sra_file, name_of_fastq_file): #For the -M flag of the fastq-dump
    '''Generates a fastq file with 10000 spots and uses this data to identify the largest read length of the dataset.
    This function is necessary to generate a full fastq with no variation in read length, a prerequisite for NOVOPlasty usage'''
    os.system("fastq-dump -X 10000 --split-spot --defline-seq '@$ac-$sn/$ri' --defline-qual '+' -O ./ %s" % (name_of_sra_file))
    with open(name_of_fastq_file) as fastq:
        length = 0
        for line in fastq:
            if line.startswith("@"):
                continue
            if len(line) > length:
                length = len(line)
            next(fastq)
            next(fastq)
        return length

def run_NOVOPlasty(accession,species,name_of_fastq_file,name_of_config_file,name_of_seed_file,max_read_length):
    with open(name_of_config_file , "a") as config: ##First, we have to prepare the configuration file.
        config.write("""Project:
-----------------------
Project name          = %s-%s
Type                  = mito
Genome Range          = 12000-22000
K-mer                 = %d
Max memory            = %d
Extended log          = 0
Save assembled reads  = no
Seed Input            = %s
Reference sequence    = 
Variance detection    = no
Heteroplasmy          = 
HP exclude list       =
Chloroplast sequence  = 

Dataset 1:
-----------------------
Read Length           = %d
Insert size           =
Platform              = illumina
Single/Paired         = PE
Combined reads        = %s
Forward reads         = 
Reverse reads         =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
""" % (species, accession, args.kmer, args.maxmemory, name_of_seed_file, max_read_length, name_of_fastq_file))
    print("Running NOVOPlasty...")
    try: ##Then, run the program(needs to be in $PATH) - Only 2.7.2 works at the moment
        os.system("NOVOPlasty2.7.2.pl -c {0}  >{1}.output 2>{1}.error &".format(name_of_config_file, name_of_config_file[:-7]))
        print("NOVOPlasty assembly successfully finished!\n")
    except:
        print("NOVOPlasty assembly error. Please check the '%s.error' file to identify the problem.\n" % (name_of_config_file[:-7]))
#Have to find a way for the command line to work with other NOVOPlasty versions (not only 2.7.2).
#Need to improve error catching (try and except). The except could be addressed by checking if anything has been written to the error file. The errors during NOVOPlasty could be caught by using tail "-n1" on the output file or by checking if the fasta sequence files have been generated.


if args.filename:
    print(main_function(args.filename))
