#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess
import os
import re
from Bio import SeqIO, Entrez
import shutil
import gzip
import traceback
from mitoattributes import mitofree_attributes

import functools
print = functools.partial(print, flush=True)

class InvalidSeedError(Exception):
    pass

class mitoassembly(mitofree_attributes): #INHERITANCE!!!
    '''This class takes the sra_run_number, species name (or unique identifier for the sample) and seed accession and performs the mitogenome assembly'''
    def __init__(self, dataset_line): #Specific assembly attributes from mitofree_attributes have to be transferred here!!!
        super().__init__(dataset_line) #When passing 'line' as argument, access the 'init' method of the base class with said line as argument. This way, all attributes from base class are available to child class and custom attributes can be added to child class too.
        
    def main_function(self):
        self.create_directory()
        if self.download_sra_files_prefetch():
            self.max_read_length = self.max_read_length()
            self.generate_fastq()
            self.download_seed()
            self.create_NOVOPlasty_config()
            self.run_NOVOPlasty()
            self.novop_result = self.check_novop_assembly()
            assert self.novop_result, "NOVOPlasty result files could not be found. Check 'novop.out' and 'novop.err'."
            self.merge_NOVOPlasty_contigs()
            if not self.check_mitobim_assembly_finished():
                self.run_mitobim()
                self.iteration = self.get_complete_iteration()
                assert self.iteration, "There must be at least one complete assembly iteration."
                self.ace = self.mitobim_convert_maf_to_ace()
                self.mitobim_ace_to_fasta()
            if self.savespace:
                self.savespace_func()

    def create_directory(self):
        if not os.path.isdir(self.prefix):
            os.mkdir(self.prefix)
        os.chdir(self.prefix)
        
    def download_sra_files_prefetch(self):
        print("WORKING DIR: {}\n".format(os.getcwd()))
        if os.path.isfile(self.sra_file): ##Checks if the sra file has already been downloaded
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
        if os.path.isfile(self.seed_file): ##Checks if the seed file has already been downloaded. NOT WORKING!
            print("This seed has already been downloaded.")
            return
        try:
            with open(self.seed_file, "w+") as fasta:
                handle = Entrez.efetch(db='nucleotide', id=self.seed, rettype='fasta', retmode='text')
                fasta.write(handle.read())
            print("Seed file %s downloaded succesfully" % (self.seed_file))
        except:
            print("Seed file %s could not be downloaded" % (self.seed_file))

    def check_file_exists(self, filename):
        if os.path.isfile(filename): ##Checks if file is available
            return True

    def generate_fastq(self):
        if self.check_file_exists(self.fastq_file) and os.stat(self.fastq_file).st_size > 100000000: #FASTQ FILE HAS TO BE LARGER THAN 100 MB. REMOVE THIS LATER.
            print("The file %s has already been converted. Assembly will proceed normally.\n" %(self.fastq_file))
            return
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
        if self.check_file_exists(self.fastq_file): #If the file is already there, it has already been normalized by length.
            with open(self.fastq_file) as fastq: #Thus, we only need to count the length of the first read (2nd line of file)
                fastq.readline() #Skips first line (header)
                seq = fastq.readline() # gets 2nd line (sequence)
                length = len(seq) 
                #print(length, seq)
            return length
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

    def check_novop_assembly(self):
        for i in [self.novop_assembly_circular, self.novop_assembly_merged, self.novop_assembly_partial]:
            if self.check_file_exists(i) and os.stat(i).st_size != 0:
                return i
        else:
            return False
    
    def create_NOVOPlasty_config(self):
        config_path = "{}/novop_config_template.txt".format(self.scriptdir)
        with open(config_path) as template, open(self.config_file, "w") as config_out:
            config_out.write(template.read().format(self.prefix, self.novop_kmer, self.maxmemory, self.seed_file, self.max_read_length, self.fastq_file))
            
    def run_NOVOPlasty(self):
        if self.check_novop_assembly():
            print("NOVOPlasty assembly already finished. Going forward...")
            return
        print("Running NOVOPlasty...")
        try:
            with open("novop.out", "w+") as output, open('novop.err', 'w+') as error:
                novop_process = subprocess.run(["NOVOPlasty3.7.2.pl", "-c", self.config_file], timeout=self.timeout*3600, stdout=output, stderr=error, check=True)
                output.seek(0)
                for line in output:
                    if line.startswith("INVALID SEED"):
                        raise InvalidSeedError("INVALID SEED: Please try a different seed sequence")
        except subprocess.TimeoutExpired:
            print("NOVOPlasty assembly taking too long. Skipping to next dataset if there is any")
            return False
        except subprocess.CalledProcessError:
            print("Something went wrong with the NOVOPlasty assembly. Check 'novop.err' for more info")
            return False
        else:
            print("NOVOPlasty assembly finished!")

    def count_contigs(self):
        number_of_contigs = 0
        print(self.novop_result)
        with open(self.novop_result) as fa:
            for line in fa:
                if line.startswith(">"):
                    number_of_contigs += 1
        return number_of_contigs
        
    def merge_NOVOPlasty_contigs(self):
        '''Uses CAP3 to merge contigs assembled by NOVOPlasty'''
        assert self.novop_result, "NOVOPlasty result files could not be found. Please check 'novop.out' and 'novop.err' for assembly errors."
        number_of_contigs = self.count_contigs()
        if number_of_contigs > 1:
            print("This assembly has %d contigs. Attempting to merge them with CAP3..." % (number_of_contigs))
            with open("cap3.out", "w") as output, open("cap3.err", "w") as error:
                cap3_alignment = subprocess.Popen(["cap3", self.novop_result], stdout=output, stderr=error)
                cap3_alignment.wait()
                self.get_largest_contig()
                ##In this section, I could add an if statement to check if the 'singlets' file is empty AND the 'contigs' file has a single contig. If it does, I should tell the user about that. Information about other cases could be added too.
        else:
            print("This assembly has only %d contig. Nothing to be merged here." % (number_of_contigs))
            with open(self.novop_result) as fasta, open("largest_contig.fa", "w") as largest:
                largest.write(fasta.read())
        return True

    def get_largest_contig(self): ##Gets the largest contig assembled by CAP3 and writes it to a file
        merged = "{}.cap.contigs".format(self.novop_result)
        unmerged = "{}.cap.singlets".format(self.novop_result)
        results = [merged, unmerged] #CAP3 outputs two files: one with the contigs assembled, and one with singlets (unassembled sequences). We need to look for the largest sequence in both file.
        max_len = 0
        full_record = ""
        for fasta in results: ##This 'for' statement gets the largest contig...
            for sequence in SeqIO.parse(fasta, "fasta"):
                if len(sequence) > max_len:
                    max_len = len(sequence)
                    full_record = ">%s\n%s" % (self.prefix, sequence.seq) #Change sequence id in order to become compatible with MITObim (Mira does not accept contigs called "contigs"
        with open("largest_contig.fa", "w") as contig:
            contig.write(full_record)

    def check_mitobim_assembly_finished(self):
        if self.check_file_exists(self.mitobim_result) and os.stat(self.mitobim_result).st_size != 0:
            print("MITObim assembly already finished. Going straight to the annotation process")
            return True
        else:
            return False            
            
    def run_mitobim(self): ##NEED TO IMPLEMENT TIMEOUT
        print("Removing any iteration folders for previous MITObim assemblies...")
        self.remove_mitobim_iterations()
        print("Running MITObim for species {}...".format(self.species))
        print("Command used: MITObim.pl -end 100 -quick largest_contig.fa -sample {} -ref mitobim -readpool {} --clean".format(self.species, self.fastq_file))
        try:
            with open("mitobim.out", "w") as out, open("mitobim.err", "w") as err:
                mitobim_path = "{}/MITObim.pl".format(self.scriptdir)
                mitobim = subprocess.run([mitobim_path, "-end", "100", "-quick", "largest_contig.fa", "-sample", self.species, "-ref", "mitobim", "-readpool", self.fastq_file, "--clean"], timeout=self.timeout*3600, stdout=out, stderr=err, check=True) ##--clean should be an optional parameter in the final version of the script
            print("MITObim assembly succesfully finished!")
            return True
        except subprocess.TimeoutExpired:
            print("MITObim assembly taking too long. Will get the last iteration's result and skip to next dataset if there is any")
            return False
        except subprocess.CalledProcessError:
            print("Something went wrong with the MITObim assembly. Check 'mitobim.err' for more info")
            return False
            
    def iteration_list(self):
        iterations = [i for i in os.listdir(".") if i.startswith("iteration")] ##Get a list of all iteration directories generated by MITObim in the current directory
        return iterations

            
    def get_mitobim_last_iteration(self):
        print(os.getcwd())
        iterations = self.iteration_list()
        print(iterations)
        max_iteration_number = 0
        for i in iterations:
            current_iteration_number = int(re.match("^iteration(\d+)$", i).group(1)) #Get the number of each iteration. If it is greater than the previous value of max_iteration_number, its value is updated and the iteration folder is saved in the "last_iteration" variable
            if current_iteration_number > max_iteration_number:
                max_iteration_number = current_iteration_number
        return(max_iteration_number) ##In the end, returns the last iteration generated by MITObim

    def get_complete_iteration(self):
        '''Checks if the final iteration has been finalized. If the last iteration has not been finalized (timeout), it works on the second last iteration directory'''
        max_iteration = self.get_mitobim_last_iteration()
        print(max_iteration)
        if max_iteration == 0:
            print("Iteration 0 has not been finished. Please look for potential issues on 'mitobim.err'. You can also try running this assembly with a subsample of the data")
            return False
        iterations = ["iteration{}".format(max_iteration), "iteration{}".format(max_iteration-1)] #Since the --clean flag is being used, only the last and second last iterations are available.
        for i in iterations:
            if os.path.isfile("{0}/{1}-mitobim_assembly/{1}-mitobim_d_results/{1}-mitobim_out.maf".format(i, self.species)):
                return(i)
        print("No iteration folder has a finalized assembly... Proceeding to the next species...")
        return False

    def gzip_ace(self):
        gz = "{}.gz".format(self.ace)
        print("Compressing ACE...")
        with open(self.ace, "rb") as consensus, gzip.open(gz, "wb") as gz:
            gz.write(consensus.read())
        print("ACE compressed to gzip. Removing original uncompressed file...")
        os.remove(self.ace)
        print("Uncompressed ACE file removed.")

    def mitobim_convert_maf_to_ace(self):
        mitobim_ref_name = "mitobim"
        mitobim_prefix = "{}-{}".format(self.species, mitobim_ref_name)
        maf = "{0}/{1}_assembly/{1}_d_results/{1}_out.maf".format(self.iteration, mitobim_prefix)
        ace = "{}.{}".format(mitobim_prefix, self.iteration).upper()
        print("Converting MAF to ACE...")
        miraconvert = subprocess.Popen(["miraconvert", "-f", "maf", "-t", "ace", "-r", "C", "-r", "f", maf, ace])
        miraconvert.wait()
        print("Conversion to ACE finished!")
        ace = "{}.ace".format(ace)
        return(ace)

    def mitobim_ace_to_fasta(self):
        print("Extracting assembly from ACE file...")
        with open(self.mitobim_result, "w") as it:
            for seq in SeqIO.parse(self.ace, "ace"):
                it.write(seq.format("fasta").replace("-", "")) #The assembly generally contains gap characters "-" that need to be removed
        print("Assembly saved to {}".format(self.mitobim_result))
        self.gzip_ace()

    def remove_mitobim_iterations(self):
        iterations = self.iteration_list()
        for i in iterations:
            shutil.rmtree(i)

    def savespace_func(self):
        os.remove(self.fastq_file)
        self.remove_mitobim_iterations()

#if args.filename:
#     print(main_function(args.filename))

'''# get rid of Biopython warning in CAI -- could it be used to remove the Entrez warning?
import warnings
from Bio import BiopythonWarning

warnings.simplefilter("ignore", BiopythonWarning)
'''
