'''
Created on Feb 8, 2012

@author: mkiyer
'''
from config import DefuseConfig, FILTERED_FASTQ_PREFIX

def get_defuse_config_string(species_dir, defuse_config,
                             max_fragment_size,
                             num_processors):
    species_sub = lambda arg: arg.replace("${SPECIES}", species_dir)
    lines = ["source_directory\t=\t%s" % (defuse_config.source_dir)]
    for k in DefuseConfig.path_params.iterkeys():
        v = getattr(defuse_config, k)
        lines.append("%s\t=\t%s" % (k,species_sub(v)))
    for k in DefuseConfig.params.iterkeys():
        v = getattr(defuse_config, k)
        lines.append("%s\t=\t%s" % (k,v))
    lines.extend(["# Paths to external tools",
                  "bowtie_bin\t=\tbowtie",
                  "bowtie_build_bin\t=\tbowtie-build",
                  "blat_bin\t=\tblat",
                  "fattotwobit_bin\t=\tfaToTwoBit",
                  "r_bin\t=\tR",
                  "rscript_bin\t=\tRscript",
                  "# Dataset files",
                  "dataset_prefix                              = $(dataset_directory)/defuse",                 
                  "chromosome_prefix                           = $(dataset_prefix).dna.chromosomes",                 
                  "exons_fasta                                 = $(dataset_prefix).exons.fa",
                  "cds_fasta                                   = $(dataset_prefix).cds.fa",
                  "cdna_regions                                = $(dataset_prefix).cdna.regions",
                  "cdna_fasta                                  = $(dataset_prefix).cdna.fa",
                  "reference_fasta                             = $(dataset_prefix).reference.fa",
                  "rrna_fasta                                  = $(dataset_prefix).rrna.fa",
                  "ig_gene_list                                = $(dataset_prefix).ig.gene.list",
                  "repeats_regions                             = $(dataset_directory)/repeats.regions",
                  "est_split_fasta1                            = $(dataset_directory)/est.1.fa",
                  "est_split_fasta2                            = $(dataset_directory)/est.2.fa",
                  "est_split_fasta3                            = $(dataset_directory)/est.3.fa",
                  "est_split_fasta4                            = $(dataset_directory)/est.4.fa",
                  "est_split_fasta5                            = $(dataset_directory)/est.5.fa",
                  "est_split_fasta6                            = $(dataset_directory)/est.6.fa",
                  "est_split_fasta7                            = $(dataset_directory)/est.7.fa",
                  "est_split_fasta8                            = $(dataset_directory)/est.8.fa",
                  "est_split_fasta9                            = $(dataset_directory)/est.9.fa",
                  "# Fasta files with bowtie indices for prefiltering reads for concordantly mapping pairs",
                  "prefilter1                                  = $(unigene_fasta)",
                  "# deFuse scripts and tools",
                  "scripts_directory                           = $(source_directory)/scripts",
                  "tools_directory                             = $(source_directory)/tools",
                  "data_directory                              = $(source_directory)/data",
                  "# Path to samtools 0.1.8 is compiled for you, use other versions at your own risk",
                  "samtools_bin                                = $(source_directory)/external/samtools-0.1.8/samtools",
                  "# Bowtie parameters",
                  "bowtie_threads                              = %d" % (num_processors),
                  "bowtie_quals                                = --phred33-quals",
                  "max_insert_size                             = %d" % (max_fragment_size),
                  "# Parameters for building the dataset",
                  "chromosomes                                 = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT",
                  "mt_chromosome                               = MT",
                  "gene_sources                                = IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,processed_transcript,protein_coding",
                  "ig_gene_sources                             = IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_pseudogene",
                  "rrna_gene_sources                           = Mt_rRNA,rRNA,rRNA_pseudogene",
                  "# Blat sequences per job",
                  "num_blat_sequences                          = 10000",
                  "# Number of reads for each job in split",
                  "reads_per_job                               = 1000000",
                  "# Number of regions for each breakpoint sequence job in split",
                  "regions_per_job                             = 20",
                  "# If you have command line 'mail' and wish to be notified",
                  "#mailto                                      = andrew.mcpherson@gmail.com",
                  "# Remove temp files",
                  "remove_job_files                            = yes",
                  "remove_job_temp_files                       = yes",
                  "# Converting to fastq",
                  "data_lane_regex_1                           = ^(%s)[12]\.fq.*$" % (FILTERED_FASTQ_PREFIX),
                  "data_end_regex_1                            = ^%s([12])\.fq.*$" % (FILTERED_FASTQ_PREFIX),
                  "data_compress_regex_1                       = ^%s[12]\.fq(.*)$" % (FILTERED_FASTQ_PREFIX),
                  "data_converter_1                            = cat"])
    return '\n'.join(lines)


class DefuseConfig(object):    
    path_params = {"gene_models": "",
                   "genome_fasta": "",
                   "repeats_filename": "",
                   "est_fasta": "",
                   "est_alignments": "",
                   "unigene_fasta": "",
                   "dataset_directory": ""}
    params = {"clustering_precision": "0.95",
              "span_count_threshold": "5",
              "split_count_threshold": "3",
              "percent_identity_threshold": "0.90",
              "max_dist_pos": "600",
              "num_dist_genes": "500",
              "split_min_anchor": "4",
              "max_concordant_ratio": "0.1",
              "splice_bias": "10",
              "denovo_assembly": "no",
              "positive_controls": "$(data_directory)/controls.txt",
              "probability_threshold": "0.50",
              "covariance_sampling_density": "0.01",
              "dna_concordant_length": "2000",
              "discord_read_trim": "50"}

    def __init__(self):
        # check environment for location of source code
        if "DEFUSEPATH" in os.environ:
            self.source_dir = os.environ["DEFUSEPATH"]
        else:
            self.source_dir = ""
        for k,v in DefuseConfig.path_params.iteritems():
            setattr(self, k, v)
        for k,v in DefuseConfig.params.iteritems():
            setattr(self, k, v)

    @staticmethod
    def from_xml_elem(elem):
        c = DefuseConfig()
        for subelem in elem.iter():
            setattr(c, subelem.tag, subelem.text)
        return c
    
    def to_xml(self, root):
        for attrname in DefuseConfig.path_params.iterkeys():
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for attrname in DefuseConfig.params.iterkeys():
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
    
    def is_valid(self, species_dir):
        valid = True
        species_sub = lambda arg: arg.replace("${SPECIES}", species_dir)
        # required paths
        for attrname in DefuseConfig.path_params.iterkeys():
            val = getattr(self, attrname)
            newval = species_sub(val)
            if not os.path.exists(newval):
                logging.error("File not found: %s" % (newval))
                valid = False
        return valid