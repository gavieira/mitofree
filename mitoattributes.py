import os

class mitofree_attributes():
    '''Base class that holds all general info for MitoFree runs'''
    def __init__(self, dataset_line, novop_kmer=39, mitob_kmer=73, maxmemory=0, subset=50000000, timeout=24, savespace=False): ##IMPLEMENT ARGS AND KWARGGS?
        self.timeout = timeout*3600 #Converting to hours
        self.scriptdir = os.path.dirname(os.path.realpath(__file__))
        self.novop_kmer = novop_kmer
        self.mitob_kmer = mitob_kmer
        self.maxmemory = maxmemory
        self.subset = subset
        self.savespace = savespace
        
        self.sra_run_number = dataset_line.split("\t")[0].strip() #SRA run number for the sequencing dataset
        self.species = dataset_line.split("\t")[1].strip()
        self.prefix = "{}-{}".format(self.species, self.sra_run_number)
        self.seed = dataset_line.split("\t")[2].strip()
        #self.new_working_dir = "%s/%s-%s" % (base_working_dir, species.upper(), accession.upper())
        self.sra_file = "{}.sra".format(self.prefix)
        self.fastq_file = "{}.fastq".format(self.prefix)
        self.config_file = "{}.config".format(self.prefix)
        self.seed_file = "{}.seed.fa".format(self.seed)
        self.novop_assembly_circular = "Circularized_assembly_1_{}-{}.fasta".format(self.species, self.sra_run_number) # The "1" should be changed to regex in order to accept any digit, but os.path.isfile (used in the "merge contigs" section) does not work with regex 
        self.novop_assembly_merged = "Option_1_{}-{}.fasta".format(self.species, self.sra_run_number) #In this case, NOVOPlasty managed to merge the contigs, and if this file contains only one contig, we are going to use it for the next steps without the use of CAP3.
        self.novop_assembly_partial = "Contigs_1_{}-{}.fasta".format(self.species, self.sra_run_number) #Partial assemblies, unmerged
        self.mitobim_result = "{}_mitobim.fasta".format(self.prefix)