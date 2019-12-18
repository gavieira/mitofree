#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess
import os
import re
import fileinput
from Bio import SeqIO, Entrez, BiopythonWarning
import functools
import shutil
import gzip
import signal
import traceback
import warnings



class mitoassembly():
    '''This class takes the sra_run_number, species name (or unique identifier for the sample) and seed accession and performs the mitogenome assembly''' 
    
#     Add *args and **kwargs to this class
#     Put create_folder function on the main script
    
    def __init__(self, dataset_line, kmer=39, maxmemory=0, subset=50000000): ##IMPLEMENT ARGS AND KWARGGS?
        self.scriptdir = os.path.dirname(os.path.realpath(__file__))
        self.kmer = kmer
        self.maxmemory = maxmemory
        self.subset = subset
        self.sra_run_number = dataset_line.split("\t")[0].strip() #SRA run number for the sequencing dataset
        self.species = dataset_line.split("\t")[1].strip()
        self.prefix = "{}-{}".format(self.species, self.sra_run_number)
        self.seed = dataset_line.split("\t")[2].strip()
        #self.new_working_dir = "%s/%s-%s" % (base_working_dir, species.upper(), accession.upper())
        self.sra_file = "{}.sra".format(self.prefix)
        self.fastq_file = "{}.fastq".format(self.prefix)
        self.config_file = "{}.config".format(self.prefix)
        self.seed_file = "{}.seed.fa".format(self.seed)
        self.novop_assembly_circular = "Circularized_assembly_1_{}-{}.fasta".format(self.species, self.sra_run_number) # The "1" should be changed to regex in order to accept any digit, but os.path.isfile (used in the "merge contigs" section) does not work with regex 
        self.novop_assembly_merged = "Option_1_{}-{}.fasta".format(self.species, self.sra_run_number) #In this case, NOVOPlasty managed to merge the contigs, and if this file contains only one contig, we are going to use it for the next steps without the use of CAP3.
        self.novop_assembly_partial = "Contigs_1_{}-{}.fasta".format(self.species, self.sra_run_number) #Partial assemblies, unmerged


    def main_function(self):
        try:
            self.create_directory()
            if self.download_sra_files_prefetch():
                self.max_read_length = self.max_read_length()
                self.generate_fastq()
                self.download_seed()
                self.create_NOVOPlasty_config()
                self.run_NOVOPlasty()
#                 if self.merge_priority():##Could use this to check if NOVOPlasty assembly has successfully finished and skip this step.
#                     print("NOVOPlasty assembly succesfully finished!")
#                     self.changeid_pre_mitobim()
#                     ##Add while loop that runs mitobim and, if timeout, run generate_fastq with half the read number and then runs mitobim again:
#                     #while not run_mitobim("largest_contig.fa", species, name_of_fastq_file):
#                     #args.subset = args.subset/2
#                     #run_mitobim("largest_contig.fa", species, name_of_fastq_file)
#                     self.run_mitobim()
#                     last_it = self.last_finalized_iteration()
#                     ace = self.mitobim_convert_maf_to_ace(species, last_it)
#                     assembly = self.mitobim_ace_to_fasta(ace)
#                     if args.savespace:
#                         self.remove_assembly_files(name_of_fastq_file)
        except Exception as error:
            os.chdir("..")
            fullerror = traceback.format_exc()
            print("An error has occurred for this assembly: {}\n\nFULL ERROR:\n\n{}\n\nProceeding to the next assembly...\n".format(error, fullerror))
            pass
        else:
            os.chdir("..")
            pass #PUT ANNOTATION MODULE HERE!!!!!


    def create_directory(self):
        if not os.path.isdir(self.prefix):
            os.mkdir(self.prefix)
        os.chdir(self.prefix)


    def download_sra_files_prefetch(self):
        print("WORKING DIR: {}\n".format(os.getcwd()))
        if os.path.isfile(self.sra_file): ##Checks if the sra file has already been downloaded.
            print("The file %s has already been downloaded. Assembly will proceed normally.\n" %(self.sra_file))
            return True
        else:
            print("Downloading {}:\n".format(self.sra_run_number))
            prefetch = subprocess.run(["prefetch", "--max-size", "900000000",  "--location", ".", "-o", self.sra_file, self.sra_run_number], capture_output=True) ##Only works with prefetch >= 2.10.0
            with open("prefetch.out", "w") as stdout:
                stdout.write(prefetch.stdout.decode())
            if prefetch.returncode == 0:
                return True
            else:
                empty_dir = os.getcwd()
                shutil.rmtree(empty_dir)
                print("The {} dataset could not be downloaded. Prefetch stderr:\n{}\n".format(self.sra_run_number, prefetch.stderr.decode()))
                return False


    def download_seed(self):
        warnings.simplefilter('ignore', BiopythonWarning)
        try:
            with open(self.seed_file, "w+") as fasta:
                handle = Entrez.efetch(db='nucleotide', id=self.seed, rettype='fasta', retmode='text')
                fasta.write(handle.read())
            print("Seed file %s downloaded succesfully" % (self.seed_file))
        except:
            print("Seed file %s could not be downloaded" % (self.seed_file))

    def generate_fastq(self):
        print("Converting {} to fastq...".format(self.sra_run_number))
        try:
            os.system("fastq-dump -M {} -X {} --split-spot --defline-seq '@$ac-$sn/$ri' --defline-qual '+' -O ./ {}".format(self.max_read_length-1,self.subset//2,self.sra_file)) #Maximum of 1 billion reads
            print("Dataset has been converted to fastq succesfully!\n")
        #fastq_name = re.sub("sra$", "fastq", sra_file)
        #with open(fastq_name, "a") as fastq:
        except:
            print("Dataset could not be converted to fastq\n")
        #use sratoolkit's fastq-dump to generate a subset(or not) of reads with header
        #use regex to extract list of read sizes "length=([0-9]+)$"
        #identify max read length and use it as parameter for fastq-dump


    def max_read_length(self): #For the -M flag of the fastq-dump
        '''Generates a fastq file with 10000 spots and uses this data to identify the largest read length of the dataset.
        This function is necessary to generate a full fastq with no variation in read length, a prerequisite for NOVOPlasty usage'''
        os.system("fastq-dump -X 10000 --split-spot --defline-seq '@$ac-$sn/$ri' --defline-qual '+' -O ./ %s" % (self.sra_file))
        with open(self.fastq_file) as fastq:
            length = 0
            for line in fastq:
                if line.startswith("@"):
                    continue
                if len(line) > length:
                    length = len(line)
                next(fastq)
                next(fastq)
            return length

    def check_NOVOPlasty_files(self):
        for i in [self.novop_assembly_circular, self.novop_assembly_merged, self.novop_assembly_partial]:
            if os.path.isfile(i) and os.stat(i).st_size != 0:
                return True
        else:
            return False
    
    def create_NOVOPlasty_config(self):
        config_path = "{}/novop_config_template.txt".format(self.scriptdir)
        with open(config_path) as template, open(self.config_file, "w") as config_out:
            config_out.write(template.read().format(self.prefix, self.kmer, self.maxmemory, self.seed_file, self.max_read_length, self.fastq_file))
            

    def run_NOVOPlasty(self):
        if self.check_NOVOPlasty_files():
            print("NOVOPlasty assembly already finished. Going forward...")
            return
        print("Running NOVOPlasty...")
        with open("novop.out", "w") as output, open('novop.err', 'w') as error:
            novop_assembly =  subprocess.Popen(["NOVOPlasty3.0.pl", "-c", self.config_file], stdout=output, stderr=error)
            novop_assembly.wait()
    #Have to find a way for the command line to work with other NOVOPlasty versions (not only 2.7.2).
    #Need to improve error catching (try and except). The except could be addressed by checking if anything has been written to the error file. The errors during NOVOPlasty could be caught by using tail "-n1" on the output file or by checking if the fasta sequence files have been generated.

    def merge_priority(self): ##Repetitive returns and statements in "except" block are not being executed. Needs to be debbuged
        mergefile = str()
        if os.path.isfile(name_of_novop_assembly_circular) and os.stat(name_of_novop_assembly_circular).st_size != 0:
            mergefile = self.novop_assembly_circular
        elif os.path.isfile(name_of_novop_assembly_merged) and os.stat(name_of_novop_assembly_merged).st_size != 0:
            mergefile = self.novop_assembly_merged
        elif os.path.isfile(name_of_novop_assembly_partial) and os.stat(name_of_novop_assembly_partial).st_size != 0:
            mergefile = self.novop_assembly_partial    
        else: ##THE SCRIPT DOES NOT EXECUTE THESE LINES OF CODE
            print("The file {}, {} or {} could not be found in the directory {}" % (self.novop_assembly_circular, self.novop_assembly_merged, self.novop_assembly_partial, os.getpwd()))
            print("NOVOPlasty assembly error. Please check the 'novop.out' and 'novop.err' files to identify the problem.\n")
            return False
        print(mergefile)
        merge_contigs(mergefile)
        return True

    def merge_contigs(self):
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

    def changeid_pre_mitobim(self): #Change sequence id in order to become compatible with MITObim (Mira does not accept contigs called "contigs"
        content = ""
        with open("largest_contig.fa", "r") as fasta: #Opens file (read-only), alters its id and saves the whole sequence to a variable
            for line in fasta:
                if line.startswith(">"):
                    line = re.sub("^>.*$", ">%s" % (self.prefix), line)
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

    def run_mitobim(self): ##NEED TO IMPLEMENT TIMEOUT
        with open("mitobim.out", "w") as output, open("mitobim.err", "w") as error:
            print("Running MITObim for species {}...".format(species))
            print("Command used: MITObim.pl -end 100 -quick largest_contig.fa -sample {} -ref mitobim -readpool {} --clean".format(largest_contig, species, name_of_fastq_file))
            time_handler = signal.signal(signal.SIGALRM, sigalrm_handler)
            signal.alarm(args.timeout*3600) # Start timer, converting hours to minutes
            try:
                mitobim = subprocess.Popen(["MITObim.pl", "-end", "100", "-quick", "largest_contig.fa", "-sample", species, "-ref", "mitobim", "-readpool", name_of_fastq_file, "--clean"], stdout=output, stderr=error) ##--clean should be an optional parameter in the final version of the script
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
        return(consensus)


    def remove_assembly_files(name_of_fastq_file):
        os.remove(name_of_fastq_file)
        iterations = [i for i in os.listdir(".") if i.startswith("iteration")] ##Repetitive code - already appears in "mitobim_last_iteration()" function. Will refactor the code in order to eliminate this (make a new function that returns this list of iteration directories).
        for i in iterations:
            shutil.rmtree(i)


# if args.filename:
#     print(main_function(args.filename))
    
'''# get rid of Biopython warning in CAI -- could it be used to remove the Entrez warning?
import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)
'''
