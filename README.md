# Mitofree (in development)

***You can use the Zenodo DOI to cite this code:***

[![DOI](https://zenodo.org/badge/171532531.svg)](https://zenodo.org/badge/latestdoi/171532531)

***Docker image available [here](https://hub.docker.com/r/gavieira/mitofree)***

---

A pipeline for automated mitochondrial genome assembly using public data.

Dependencies:

* [CAP3](http://seq.cs.iastate.edu/cap3.html)
* [NOVOPlasty3.7.2](https://github.com/ndierckx/NOVOPlasty)
* [MIRA4.0.2](https://ufpr.dl.sourceforge.net/project/mira-assembler/MIRA/stable/mira_4.0.2_linux-gnu_x86_64_static.tar.bz2)
* [MITOS](https://gitlab.com/Bernt/MITOS)
  * Python2 and other software
* [Python3](https://www.python.org/)
  * [Biopython](https://biopython.org/)
* [sratoolkit2.10.0](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)

All of the above dependencies can be easily installed through [Bioconda](https://bioconda.github.io/). However, since Mitofree needs both Python3 and Python2 to run, manual creation of a conda environment for the package can be a little tricky. Thus, we encourage the use of the [***docker image***](https://hub.docker.com/r/gavieira/mitofree) to run this software.

## Running Mitofree with docker:

#### 1 - Install [docker](https://docs.docker.com/install/)

#### 2 - Download the latest image:

```
docker pull gavieira/mitofree:latest
```

#### 3- Use the following command to generate a container:

```
docker run --name mitofree -i -t -v ~:/mnt -w /mnt gavieira/mitofree /bin/bash
```

***OBS***: After creating the container, you will not need to go through steps 1 and 2 again. You can simply start the container anytime you want. To do so, run:

```
docker start -i mitofree
```

#### 4- Finally, run Mitofree:

```
usage: mitofree.py [-h] [-S] [-M] [--novop_kmer] [--mitob_kmer] [-g] [-s] [-T]
                   FILENAME

Downloads sra NGS data and assembles mitochondrial contigs using NOVOPlasty
and MITObim

positional arguments:
  FILENAME           Path to file with multiple accessions (one per line)

optional arguments:
  -h, --help         show this help message and exit
  -S, --savespace    Automatically removes residual assembly files such as
                     fastq and mitobim iterations
  -M , --maxmemory   Limit of RAM usage for NOVOPlasty. Default: no limit
  --novop_kmer       K-mer used in NOVOPlasty assembly. Default: 39
  --mitob_kmer       K-mer used in MITObim assembly. Default: 73
  -g , --gencode     Genetic code table. Default: 2 (Vertebrate Mitochondrial)
  -s , --subset      Max number of reads used in the assembly process.
                     Default: 50 million reads
  -T , --timeout     Custom timeout for MITObim, in hours. Default: 24h
```

Please note the -M "--maxmemory" argument, that limits NOVOPlasty's RAM usage (in GB). If you are running this software from a machine with limited RAM available, you will want to set this option so that it won't use all your memory. For instance, if you have a 8GB computer, you may want to use "-M 7".

The -s "--subset" argument can be used to limit dataset size, which can also reduce RAM requirements. This argument can also be used to increase dataset size, which may be useful if you're having trouble in circularizing a mitogenome and some RAM to spare.



## Example of Mitofree's input file:

Basically, this file consists of three tab-separated collumns, each with a specific information:

1-SRA_RUN_NUMBER        2-SPECIES_NAME          3-SEED_GENBANK_ACCESSION

For instance:

```
ERR1306022	Species1	MK297287
ERR7295165	Species2	MK297241
ERR1306034	Species3	MK291745
#SRR4409513	Species4	MK291678 #This assembly will be skipped
```

Each line corresponds to a different assembly. This way, you can build a list of as many organisms as you want and assemble their mitogenomes all at once. It is also possible to skip an assembly by adding a hash symbol (#) at the start of its corresponding line.
