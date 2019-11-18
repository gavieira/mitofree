# MitoFree (in development)

***You can use the Zenodo DOI to cite this code:*** 

[![DOI](https://zenodo.org/badge/171532531.svg)](https://zenodo.org/badge/latestdoi/171532531)

---

A pipeline for automated mitochondrial genome assembly using public data.

Needs Biopython module installed for Python3. You can easily install them through pip3 or [conda](https://docs.conda.io/en/latest/).

```
#If you don't have pip3 installed,run this command (ubuntu):
sudo apt install python3-pip

#Then, you can install Biopython using:
pip3 install Biopython
```

You will also need to download and unpack [NOVOPlasty3.0](https://github.com/ndierckx/NOVOPlasty), [sratoolkit (>=2.10.0)](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) and [CAP3](http://seq.cs.iastate.edu/cap3.html), [MITObim1.9](https://github.com/chrishah/MITObim) and [MIRA4.0.2](https://ufpr.dl.sourceforge.net/project/mira-assembler/MIRA/stable/mira_4.0.2_linux-gnu_x86_64_static.tar.bz2) as well as install them in your PATH environment variable.

In order to install to the PATH variable, please add the complete path to the binaries (sometimes located in a '/bin' directory instead of the root folder for the program) of all these dependencies to your ~/.bashrc file:

```
echo 'export PATH="$PATH:/path/to/NOVOPlasty:/path/to/sratoolkit:/path/to/CAP3:/path/to/MITObim:/path/to/MIRA"' >> ~/.bashrc
```

Then, source the file (or just restart the terminal)

```
source ~/.bashrc
```

If you can use bash autocompletion to call these scripts, you have succesfully added them to the PATH variable.

***A friendly reminder:*** Most of these dependecies are available at the [bioconda](https://bioconda.github.io/) channel and thus can be easily installed through [conda](https://docs.conda.io/en/latest/). We are planning to add mitofree to bioconda in the future, which will make its instalation way simpler.


Lastly, you should download MitoFree and give it execute permission:

```
chmod +x /path/to/mitofree.py
```

And then run it:

```
/path/to/mitofree.py [-h] [-S] [-M] [-K] dataset_list.txt

optional arguments:
  -h, --help         show this help message and exit
  -S, --savespace    Automatically removes residual assembly files such as
                     fastq and mitobim iterations
  -M , --maxmemory   Limit of RAM usage for NOVOPlasty. Default: no limit
  -K , --kmer        K-mer used in NOVOPlasty assembly. Default: 39
```

Please note the -M "--maxmemory" argument, that limits NOVOPlasty's RAM usage (in GB). If you are running this software from a machine with limited RAM available, you will want to set this option so that it won't use all your memory. For instance, if you have a 8GB computer, you may want to use "-M 7". 

The 'dataset_list.txt' is a plain text file that contains three tab-separated collumns, each corresponding to a specific information used by mitofree:

1-SRA_RUN_NUMBER        2-SPECIES_NAME          3-SEED_GENBANK_ACCESSION

For instance:

```
ERR1306022	Species1	MK297287
ERR7295165	Species2	MK297241
ERR1306034	Species3	MK291745
SRR4409513	Species4	MK291678
```

Each line corresponds to a different assembly. This way, you can build a list of as many organisms as you want and assemble their mitogenomes all at once.
