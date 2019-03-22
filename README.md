# MitoFree (in development)

A pipeline for automated mitochondrial genome assembly using public data.

Needs Biopython and wget modules installed for Python3. You can easily install them through pip3.

```
#If you don't have pip3 installed,run this command (ubuntu):
sudo apt install python3-pip

#Then, you can install both wget and Biopython using:
pip3 install Biopython && pip3 install wget
```

You will also need to download and unpack [NOVOPlasty2.7.2](https://github.com/ndierckx/NOVOPlasty) and [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/), as well as install them in your PATH environment variable.

In order to install to the PATH variable, please add the complete path to NOVOPlasty and sratoolkit to your ~/.bashrc file:

```
echo 'export PATH="$PATH:/path/to/NOVOPlasty:path/to/sratoolkit/bin"' >> ~/.bashrc
```

Then, source the file (or just restart the terminal)

```
source ~/.bashrc
```

If you can use bash autocompletion to call these scripts, you have succesfully added them to the PATH variable.

Lastly, you should download MitoFree and give it execute permission:

```
chmod +x /path/to/mitofree.py
```

And then run it:

```
/path/to/mitofree.py dataset_list.txt
```

Please note the -M "--maxmemory" argument, that limits NOVOPlasty's RAM usage (in GB). If you are running this software from a machine with limited RAM available, you will want to set this option so that it won't use all your memory. For instance, if you have a 8GB computer, you may want to use "-M 7". 

The 'dataset_list.txt' is a plain text file that contains three tab-separated collumns, each corresponding to a specific information used by mitofree:

1-SRA_RUN_NUMBER  2-SPECIES_NAME  3-SEED_GENBANK_ACCESSION

For instance:

```
ERR1306022	Species1	MK297287
ERR7295165	Species2	MK297241
ERR1306034	Species3	MK291745
SRR4409513	Species4	MK291678
```

Each line corresponds to a different assembly. This way, you can build a list of as many organisms as you want and assemble their mitogenomes all at once.
