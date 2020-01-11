#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py


import argparse, functools, sys
import mitoassembly, mitoannotate

assert ('linux' in sys.platform), "This code runs on Linux only."

print = functools.partial(print, flush=True) #All "print" functions have flush=True as default. This way, its contents are not buffered, being instead flushed to the standard output. With this, stdout and stderr redirection works like a charm...


# def create_folders(new_working_dir): ##Creates folder for each dataset and changes the working directory
#         print("New working directory is '%s'\n" % (new_working_dir))
#         try:
#             output_mkdir = os.system("mkdir -p %s" % (new_working_dir))
#             os.chdir(new_working_dir)
#             return True
#         except re.search(".*Permission denied$", output_mkdir):
#             print("Could not create folder %s. Permission denied." % (new_working_dir.split("/")[-1]))
#             return False
#         except:
#             print("Could not create folder %s" % (new_working_dir.split("/")[-1]))
#             return False

# def main_function(sra_list):
#     try:
#         self.create_folders()
#         if self.download_sra_files_prefetch():
#             max_read_length = self.highest_read_length()
#             self.generate_fastq()
#             self.download_seed()
#             self.run_NOVOPlasty()
#             if self.merge_priority():##Could use this to check if NOVOPlasty assembly has successfully finished and skip this step.
#                 print("NOVOPlasty assembly succesfully finished!")
#                 self.changeid_pre_mitobim()
#                 ##Add while loop that runs mitobim and, if timeout, run generate_fastq with half the read number and then runs mitobim again:
#                 #while not run_mitobim("largest_contig.fa", species, name_of_fastq_file):
#                 #args.subset = args.subset/2
#                 #run_mitobim("largest_contig.fa", species, name_of_fastq_file)
#                 self.run_mitobim()
#                 last_it = self.last_finalized_iteration()
#                 ace = self.mitobim_convert_maf_to_ace(species, last_it)
#                 assembly = self.mitobim_ace_to_fasta(ace)
#                 if args.savespace:
#                     self.remove_assembly_files(name_of_fastq_file)
#     except Exception as error:
#         print("\nAn error has occurred for this assembly:\n\n{}\n\nProceeding to the next assembly\n".format(error))
#         continue
#     else:
#         pass #PUT ANNOTATION MODULE HERE!!!!!

# return("All done!")
# ##Add the merge contigs and count contigs here (with its ifs, for readability)        


def getArgs():
    parser = argparse.ArgumentParser(description="Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty and MITObim")
    parser.add_argument("-S", "--savespace", action="store_true", default=False, help="Automatically removes residual assembly files such as fastq and mitobim iterations")
    parser.add_argument("-M", "--maxmemory", type=int, metavar="", default=0, help="Limit of RAM usage for NOVOPlasty. Default: no limit")
    parser.add_argument("--novop_kmer", type=int, metavar="", default=39, help="K-mer used in NOVOPlasty assembly. Default: 39")
    parser.add_argument("--mitob_kmer", type=int, metavar="", default=73, help="K-mer used in MITObim assembly. Default: 73")
    parser.add_argument("-s", "--subset", type=int, metavar="", default=50000000, help="Max number of reads used in the assembly process. Default: 50 million reads")
    #parser.add_argument("-P", "--parallel", type=int, metavar="", default=1, help="Number of parallel assemblies. Default: 1")
    parser.add_argument("-T", "--timeout", type=int, metavar="", default=24, help="Custom timeout for MITObim, in hours. Default: 24h")
    parser.add_argument("filename", type=str, metavar="FILENAME", help="Path to file with multiple accessions (one per line)")
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    with open(args.filename) as datasets:
        for line in datasets:
            if line.startswith('#'): pass
            else:
                assembly = mitoassembly.mitoassembly(line)
                assembly.main_function()
                annotation = mitoannotation.mitoannotation(assembly.mitobim_result)

