#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Converts SRA Runinfo Tables (downloadable from SRA Run Selector) into MitoFree tables. Note that seed accessions still need to be manually added to the output.")
parser.add_argument("table", type=str, metavar="table", help="Path to Runinfo table (must be a comma-separated values table)")
args = parser.parse_args()

df = pd.read_csv(args.table)
for run, organism in zip(df.Run, df.Organism):
    print("{}\t{}".format(run, organism.replace(" ", "_")))
                        
