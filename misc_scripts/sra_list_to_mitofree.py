#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Converts SRA Runinfo Tables (downloadable from SRA Run Selector) into MitoFree tables. Note that seed accessions still need to be manually added to the output.")
parser.add_argument("-t", "--tab", action="store_true", default=False, help="Converts tab-delimited tables (Default: comma-separated tables)")
parser.add_argument("table", type=str, metavar="table", help="Path to Runinfo table")
args = parser.parse_args()

if args.tab:
    df = pd.read_csv(args.table, sep='\t')
else:
    df = pd.read_csv(args.table)

zipped = zip(df.Run, df.Organism)


for run, organism in sorted(zipped, key = lambda x: x[1]):
    print("{}\t{}".format(run, organism.replace(" ", "_")))
                        
