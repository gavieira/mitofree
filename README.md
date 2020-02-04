# MitoFree (in development)

***You can use the Zenodo DOI to cite this code:*** 

[![DOI](https://zenodo.org/badge/171532531.svg)](https://zenodo.org/badge/latestdoi/171532531)

***Docker image available [here](https://hub.docker.com/repository/docker/gavieira/mitofree/general)***

---

A pipeline for automated mitochondrial genome assembly using public data.

Needs Biopython module installed for Python3. You can easily install it through pip3 or [conda](https://docs.conda.io/en/latest/).


```
#To install through conda:
conda install mitos=2.0.3 biopython=1.73 mira=4.0.2 sra-tools=2.10.0 cap3=10.2011 NOVOPlasty=3.7.2 -c bioconda -m -n mitofree

#Then, you need to activate mitofree's environment:
conda activate mitofree

#OBS: You will need python>=3.7 to run this script (because, among other things, of the 'capture_output' option from subprocess.run)
#If python version is not compatible, you can update it with:
conda update python
```

You will also need to download and unpack [NOVOPlasty3.0](https://github.com/ndierckx/NOVOPlasty), [sratoolkit (>=2.10.0)](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) and [CAP3](http://seq.cs.iastate.edu/cap3.html) and [MIRA4.0.2](https://ufpr.dl.sourceforge.net/project/mira-assembler/MIRA/stable/mira_4.0.2_linux-gnu_x86_64_static.tar.bz2) as well as install them in your PATH environment variable.

In order to install to the PATH variable, please add the complete path to the binaries (sometimes located in a '/bin' directory instead of the root folder for the program) of all these dependencies to your ~/.bashrc file:

```
echo 'export PATH="$PATH:/path/to/NOVOPlasty:/path/to/sratoolkit:/path/to/CAP3:/path/to/MITObim:/path/to/MIRA"' >> ~/.bashrc
```

Then, source the file (or just restart the terminal)

```
source ~/.bashrc
```

If you can use bash autocompletion to call these scripts, you have succesfully added them to the PATH variable.

***A friendly reminder:*** Most of these dependecies are available at the [bioconda](https://bioconda.github.io/) channel and thus can be easily installed through [conda](https://docs.conda.io/en/latest/). We are planning to add mitofree to bioconda in the future, which will make its instalation way simpler. For now, it is possible to use the [docker image](https://hub.docker.com/repository/docker/gavieira/mitofree/general) or to create an mitofree environment to easily run the application in a container, without the need to install any dependencies.


Lastly, you should download MitoFree and give it execute permission:

```
chmod +x /path/to/mitofree.py
```

## Running from docker image (needs docker software installed):

Download the latest image:

```
docker pull gavieira/mitofree:latest
```

Then, go to the directory where the MitoFree's input file is and use:

```
docker run --name mitofree -i -t -v $PWD:/mnt -w /mnt gavieira/mitofree /bin/bash
```

And then run it:

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


## Example of MitoFree's input file:

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
