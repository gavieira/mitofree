#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py


import argparse, functools, sys, os
import traceback
from mitoassembly import mitoassembly
from mitoannotate import mitoannotation

assert ('linux' in sys.platform), "This code runs on Linux only."

print = functools.partial(print, flush=True) #All "print" functions have flush=True as default. This way, its contents are not buffered, being instead flushed to the standard output. With this, stdout and stderr redirection works like a charm...



def getArgs():
    parser = argparse.ArgumentParser(description="Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty and MITObim")
    parser.add_argument("--mincontigsize", type=int, metavar="", default=0, help="Minimum mitogenome contig size allowed to go through annotation process. Default: 10 kbp")
    parser.add_argument("-S", "--savespace", action="store_true", default=False, help="Automatically removes residual assembly files such as fastq and mitobim iterations")
    parser.add_argument("-M", "--maxmemory", type=int, metavar="", default=0, help="Limit of RAM usage for NOVOPlasty. Default: no limit")
    parser.add_argument("--novop_kmer", type=int, metavar="", default=39, help="K-mer used in NOVOPlasty assembly. Default: 39")
    parser.add_argument("--mitob_kmer", type=int, metavar="", default=73, help="K-mer used in MITObim assembly. Default: 73")
    parser.add_argument("-g", "--gencode", type=int, metavar="", default=2, help="Genetic code table. Default: 2 (Vertebrate Mitochondrial)")
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
                try:
                    assembly = mitoassembly(line)
                    assembly.main_function()
                    annotation = mitoannotation(line, gencode=5)
                    if annotation.contigsize >= args.mincontigsize:
                        print("Contig size: {} Proceeding to annotation...".format(annotation.contigsize))
                        annotation.run_mitos()
                        annotation.generate_gbk()
                        annotation.copy_gbk_to_new_directory()
                    os.chdir("..")
                except Exception as error:
                    os.chdir("..")
                    fullerror = traceback.format_exc()
                    print("An error has occurred for this assembly: {}\n\nFULL ERROR:\n\n{}\n\nProceeding to the next assembly...\n".format(error, fullerror))
                    pass

