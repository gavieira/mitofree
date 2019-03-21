# MitoFree (in development)

A pipeline for automated mitochondrial genome assembly using public data.

Needs Biopython and wget modules installed for Python3. You can easily install them through pip3.

```
#If you don't have pip3 installed,run this command (ubuntu):
sudo apt install python3-pip

#Then, you can install both wget and Biopython using:
pip3 install Biopython && pip3 install wget
```

You will also need [NOVOPlasty2.7.2](https://github.com/ndierckx/NOVOPlasty) and [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/) installed in your PATH environment variable.

In order to do this, please add the complete path to NOVOPlasty and sratoolkit to your ~/.bashrc file:

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

An example file for 'dataset_list.txt' is available in the package. Build list for your organism(s) of interest based on it.
