FROM continuumio/miniconda3:latest
MAINTAINER gavieira <gabrieldeusdeth@gmail.com>

##Cloning MitoFree repo into root directory and installing biopython

RUN conda install git biopython && \
    git clone https://github.com/gavieira/mitofree && \
    conda remove git

##Adding bioconda to conda software channels

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

##Creating a environment with most MitoFree dependencies and setting it as default env

RUN  conda install mitos=2.0.3 mira=4.0.2 sra-tools=2.10.0 cap3=10.2011 NOVOPlasty=3.7.2 \
     -c bioconda -m -n mitofree 
RUN echo "source activate mitofree" > ~/.bashrc

##Adding python3 from base to mitofree env 
##(MITOS annotation tool requires python2, hence this additional step)

RUN ln -s /opt/conda/bin/python /opt/conda/envs/mitofree/bin/python3

##Adding mitofree repo and env to PATH variable and setting LC_ALL (needed for mira)

ENV PATH /mitofree/:/opt/conda/envs/mitofree/bin:$PATH
ENV LC_ALL C

##Finally, containers will open up at the /mnt directory, where volumes should be mounted

WORKDIR /mnt
