#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"


import argparse ##Put argparse in function and use __name__ == "__main__"
parser = argparse.ArgumentParser(description="Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty and MITObim")
parser.add_argument("-S", "--savespace", action="store_true", default=False, help="Automatically removes residual assembly files such as fastq and mitobim iterations")
parser.add_argument("-M", "--maxmemory", type=int, metavar="", default=0, help="Limit of RAM usage for NOVOPlasty. Default: no limit")
parser.add_argument("-K", "--kmer", type=int, metavar="", default=39, help="K-mer used in NOVOPlasty assembly. Default: 39")
parser.add_argument("-s", "--subset", type=int, metavar="", default=50000000, help="Max number of reads used in the assembly process. Default: 50 million reads")
#parser.add_argument("-P", "--parallel", type=int, metavar="", default=1, help="Number of parallel assemblies. Default: 1")
parser.add_argument("-T", "--timeout", type=int, metavar="", default=24, help="Custom timeout for MITObim, in hours. Default: 24h")
parser.add_argument("filename", type=str, metavar="FILENAME", help="Path to file with multiple accessions (one per line)")
args = parser.parse_args()

import subprocess
import os
import re
import fileinput
from Bio import SeqIO, Entrez
import functools
import shutil
import gzip
import signal

print = functools.partial(print, flush=True) #All "print" functions have flush=True as default. This way, its contents are not buffered, being instead flushed to the standard output. With this, stdout and stderr redirection works like a charm...

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
            try:
                if create_folders(new_working_dir):
                    if download_sra_files_prefetch(accession, name_of_sra_file, new_working_dir):
                        max_read_length = highest_read_length(name_of_sra_file, name_of_fastq_file)
                        generate_fastq(name_of_sra_file, max_read_length)
                        download_seed(name_of_seed_file, seed)
                        run_NOVOPlasty(accession, species, name_of_fastq_file, name_of_config_file, name_of_seed_file,max_read_length, name_of_novop_assembly_circular, name_of_novop_assembly_merged, name_of_novop_assembly_partial)
                        if merge_priority(name_of_novop_assembly_circular, name_of_novop_assembly_merged, name_of_novop_assembly_partial, new_working_dir):##Could use this to check if NOVOPlasty assembly has successfully finished and skip this step.
                            print("NOVOPlasty assembly succesfully finished!")
                            changeid_pre_mitobim("largest_contig.fa", "{}-{}".format(species, accession))
                            ##Add while loop that runs mitobim and, if timeout, run generate_fastq with half the read number and then runs mitobim again:
                            #while not run_mitobim("largest_contig.fa", species, name_of_fastq_file):
                            #args.subset = args.subset/2
                            #run_mitobim("largest_contig.fa", species, name_of_fastq_file)
                            run_mitobim("largest_contig.fa", species, name_of_fastq_file)
                            last_it = last_finalized_iteration(species)
                            ace = mitobim_convert_maf_to_ace(species, last_it)
                            mitobim_ace_to_fasta(ace)
                            if args.savespace:
                                remove_assembly_files(name_of_fastq_file)
            except Exception as error:
                print("\nAn error has occurred for this assembly:\n\n{}\n\nProceeding to the next assembly\n".format(error))
                continue
            else:
                pass #PUT ANNOTATION MODULE HERE!!!!!
        return("All done!")
##Add the merge contigs and count contigs here (with its ifs, for readability)

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

def download_sra_files_prefetch(accession, name_of_sra_file, new_working_dir):
    if os.path.isfile(name_of_sra_file): ##Checks if the sra file has already been downloaded.
        print("The file %s has already been downloaded. Assembly will proceed normally.\n" %(name_of_sra_file))
        return True
    else:
        try:
            print("Downloading %s:" % (accession))
            prefetch = subprocess.Popen(["prefetch", "--max-size", "900000000",  "--location", ".", "-o", name_of_sra_file, accession]) ##Only works with prefetch >= 2.10.0
            prefetch.wait()
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
        os.system("fastq-dump -M {} -X {} --split-spot --defline-seq '@$ac-$sn/$ri' --defline-qual '+' -O ./ {}".format(max_read_length-1,args.subset//2,name_of_sra_file)) #Maximum of 1 billion reads
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

def check_NOVOP_files(name_of_novop_assembly_circular, name_of_novop_assembly_merged,name_of_novop_assembly_partial):
    for i in [name_of_novop_assembly_circular,name_of_novop_assembly_merged, name_of_novop_assembly_partial]:
        if os.path.isfile(i) and os.stat(i).st_size != 0:
            return True
    else:
        return False
    
def run_NOVOPlasty(accession, species, name_of_fastq_file, name_of_config_file, name_of_seed_file, max_read_length, name_of_novop_assembly_circular,  name_of_novop_assembly_merged, name_of_novop_assembly_partial): ##Put the config as a separate file. Add a function to just read and modify its contents.
    if check_NOVOP_files(name_of_novop_assembly_circular, name_of_novop_assembly_merged,name_of_novop_assembly_partial):
        print("NOVOPlasty assembly already finished. Going forward...")
        return
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

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
""" % (species, accession, args.kmer, args.maxmemory, name_of_seed_file, max_read_length, name_of_fastq_file))
    print("Running NOVOPlasty...")
    with open("novop.out", "w") as output, open('novop.err', 'w') as error:
        novop_assembly =  subprocess.Popen(["NOVOPlasty3.0.pl", "-c", name_of_config_file], stdout=output, stderr=error)
        novop_assembly.wait()
#Have to find a way for the command line to work with other NOVOPlasty versions (not only 2.7.2).
#Need to improve error catching (try and except). The except could be addressed by checking if anything has been written to the error file. The errors during NOVOPlasty could be caught by using tail "-n1" on the output file or by checking if the fasta sequence files have been generated.

def merge_priority(name_of_novop_assembly_circular, name_of_novop_assembly_merged, name_of_novop_assembly_partial, new_working_dir): ##Repetitive returns and statements in "except" block are not being executed. Needs to be debbuged
    mergefile = str()
    if os.path.isfile(name_of_novop_assembly_circular) and os.stat(name_of_novop_assembly_circular).st_size != 0:
        mergefile = name_of_novop_assembly_circular
    elif os.path.isfile(name_of_novop_assembly_merged) and os.stat(name_of_novop_assembly_merged).st_size != 0:
        mergefile = name_of_novop_assembly_merged
    elif os.path.isfile(name_of_novop_assembly_partial) and os.stat(name_of_novop_assembly_partial).st_size != 0:
        mergefile = name_of_novop_assembly_partial    
    else: ##THE SCRIPT DOES NOT EXECUTE THESE LINES OF CODE
        print("The file %s, %s or %s could not be found in the directory (new_working_dir is %s)" % (name_of_novop_assembly_circular, name_of_novop_assembly_merged, name_of_novop_assembly_partial, new_working_dir))
        print("NOVOPlasty assembly error. Please check the 'novop.out' and 'novop.err' files to identify the problem.\n")
        return False
    print(mergefile)
    merge_contigs(mergefile)
    return True

def merge_contigs(name_of_novop_assembly):
    '''Uses CAP3 to merge contigs assembled by NOVOPlasty'''
    number_of_contigs = count_contigs(name_of_novop_assembly)
    if number_of_contigs > 1: ##
        print("This assembly has %d contigs. Attempting to merge them with CAP3..." % (number_of_contigs))
        with open("cap3.out", "w") as output, open("cap3.err", "w") as error:
            cap3_alignment = subprocess.Popen(["cap3", name_of_novop_assembly], stdout=output, stderr=error)
            cap3_alignment.wait()
            ##In this section, I could add an if statement to check if the 'singlets' file is empty AND the 'contigs' file has a single contig. If it does, I should tell the user about that. Information about other cases could be added too.
            get_largest_contig(name_of_novop_assembly)
    else:
        print("This assembly has only %d contig. Nothing to be merged here." % (number_of_contigs))
        with open(name_of_novop_assembly) as fasta, open("largest_contig.fa", "w") as largest:
            largest.write(fasta.read())

def count_contigs(fasta):
    with open(fasta) as fa:
        number_of_contigs = 0
        for line in fa:
            if line.startswith(">"):
                number_of_contigs += 1
    return number_of_contigs

def get_largest_contig(name_of_novop_assembly): ##Gets the largest contig assembled by CAP3 and writes it to a file
    merged = "{}.cap.contigs".format(name_of_novop_assembly)
    unmerged = "{}.cap.singlets".format(name_of_novop_assembly)
    results = [merged, unmerged] #CAP3 outputs two files: one with the contigs assembled, and one with singlets (unassembled sequences). We need to look for the largest sequence in both file.
    max_len = 0
    full_record = ""
    for fasta in results: ##This 'for' statement gets the largest contig...
        for sequence in SeqIO.parse(fasta, "fasta"):
            if len(sequence) > max_len:
                max_len = len(sequence)
                full_record = ">%s\n%s" % (sequence.id, sequence.seq)
    with open("largest_contig.fa", "w") as contig:
        contig.write(full_record)

def changeid_pre_mitobim(largest_contig, species_and_accession): #Change sequence id in order to become compatible with MITObim (Mira does not accept contigs called "contigs"
    content = ""
    with open(largest_contig, "r") as fasta: #Opens file (read-only), alters its id and saves the whole sequence to a variable
        for line in fasta:
            if line.startswith(">"):
                line = re.sub("^>.*$", ">%s" % (species_and_accession), line)
                content += line
            else:
                content += line
    with open(largest_contig, "w") as fasta: #Opens the same file (write mode), overwriting it with the contents of the variable
        fasta.write(content)
        
class TimeoutException(Exception): # Custom exception for the timeout
    pass

def sigalrm_handler(signum, frame):# Handler function to be called when SIGALRM is received
    # If we get signal, a TimeoutException will be raised.
    raise TimeoutException()
        
def run_mitobim(largest_contig, species, name_of_fastq_file): ##NEED TO IMPLEMENT TIMEOUT
    with open("mitobim.out", "w") as output, open("mitobim.err", "w") as error:
        print("Running MITObim for species {}...".format(species))
        print("Command used: MITObim.pl -end 100 -quick {} -sample {} -ref mitobim -readpool {} --clean".format(largest_contig, species, name_of_fastq_file))
        time_handler = signal.signal(signal.SIGALRM, sigalrm_handler)
        signal.alarm(args.timeout*3600) # Start timer, converting hours to minutes
        try:
            mitobim = subprocess.Popen(["MITObim.pl", "-end", "100", "-quick", largest_contig, "-sample", species, "-ref", "mitobim", "-readpool", name_of_fastq_file, "--clean"], stdout=output, stderr=error) ##--clean should be an optional parameter in the final version of the script
            mitobim.wait()
            print("MITObim assembly succesfully finished!")
            return True
        except TimeoutException:
            print("MITObim assembly taking too long. Will get the last iteration's result and skip to next dataset")
            return False
        finally:
            signal.alarm(0) # Turn off timer
            signal.signal(signal.SIGALRM, time_handler) # Restore handler to previous value

def mitobim_last_iteration():
    iterations = [i for i in os.listdir(".") if i.startswith("iteration")] ##Get a list of all iteration directories generated by MITObim in the current directory
    max_iteration_number = 0
    for i in iterations:
        current_iteration_number = int(re.match("^iteration(\d+)$", i).group(1)) #Get the number of each iteration. If it is greater than the previous value of max_iteration_number, its value is updated and the iteration folder is saved in the "last_iteration" variable
        if current_iteration_number > max_iteration_number:
            max_iteration_number = current_iteration_number
    return(max_iteration_number) ##In the end, returns the last iteration generated by MITObim

def last_finalized_iteration(species):
    '''Checks if the final iteration has been finalized. If the last iteration has not been finalized (timeout), it works on the second last iteration directory'''
    last_it = mitobim_last_iteration()
    if last_it == 0:
        print("Iteration 0 has not been finished. Please look for potential issues on '{}/mitofree.err'. You can also try running this assembly with a subsample of the data")
        return False
    iterations = ["iteration{}".format(last_it), "iteration{}".format(last_it-1)] #Since the --clean flag is being used, only the last and second last iterations are available.
    for i in iterations:
        if os.path.isfile("{0}/{1}-mitobim_assembly/{1}-mitobim_d_results/{1}-mitobim_out.maf".format(i, species)):
            return(i)
    print("No iteration folder has a finalized assembly... Proceeding to the next species...")
    return False
            
def gzip_ace(ace):
    gz = "{}.gz".format(ace)
    print("Compressing ACE...")
    with open(ace, "rb") as consensus, gzip.open(gz, "wb") as gz:
        gz.write(consensus.read())
    print("ACE compressed to gzip. Removing original uncompressed file...")
    os.remove(ace)
    print("Uncompressed ACE file removed.")
    
def mitobim_convert_maf_to_ace(species, iteration):
    ref = "mitobim"
    mitobim_prefix = "{}-{}".format(species, ref)
    maf = "{0}/{1}_assembly/{1}_d_results/{1}_out.maf".format(iteration, mitobim_prefix)
    ace = "{}.{}".format(mitobim_prefix, iteration).upper()
    print("Converting MAF to ACE...")
    miraconvert = subprocess.Popen(["miraconvert", "-f", "maf", "-t", "ace", "-r", "C", "-r", "f", maf, ace])
    miraconvert.wait()
    print("Conversion to ACE finished!")
    ace = "{}.ace".format(ace)
    return(ace)

def mitobim_ace_to_fasta(ace):
    consensus = os.path.splitext(ace)[0] + ".fasta"
    print("Extracting assembly from ACE file...")
    with open(consensus, "w") as it:
        for seq in SeqIO.parse(ace, "ace"):
            it.write(seq.format("fasta").replace("-", "")) #The assembly generally contains gap characters "-" that need to be removed
    print("Assembly saved to {}".format(consensus))
    gzip_ace(ace)


def remove_assembly_files(name_of_fastq_file):
    os.remove(name_of_fastq_file)
    iterations = [i for i in os.listdir(".") if i.startswith("iteration")] ##Repetitive code - already appears in "mitobim_last_iteration()" function. Will refactor the code in order to eliminate this (make a new function that returns this list of iteration directories).
    for i in iterations:
        shutil.rmtree(i)


if args.filename:
    print(main_function(args.filename))
    
'''# get rid of Biopython warning in CAI -- could it be used to remove the Entrez warning?
import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)
'''
