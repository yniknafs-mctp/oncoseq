'''
Created on Aug 3, 2011

@author: mkiyer
'''
import os
import logging
import xml.etree.cElementTree as etree

from base import check_executable, indent_xml
from seqdb import Patient, \
    SAMPLE_TYPE_EXOME_TUMOR, SAMPLE_TYPE_EXOME_NORMAL, \
    SAMPLE_TYPE_RNASEQ, SAMPLE_TYPE_CAPTURE_RNASEQ

from oncoseq.rnaseq.lib import validate as validate_rnaseq 

# job return codes
JOB_SUCCESS = 0
JOB_ERROR = 1

# XML files
PATIENT_XML_FILE = "patient.xml"
SAMPLE_GROUP_XML_FILE = "sample_group.xml"
SAMPLE_XML_FILE = "sample.xml"

# remote job constants
REMOTE_ANALYSIS_XML_FILE = "analysis.xml"
REMOTE_CONFIG_XML_FILE = "config.xml"
REMOTE_CODE_TARGZ_FILE = "code.tar.gz"
REMOTE_COPY_MAX_SIZE_BYTES = (8 << 30)

# fastq file names
READ1_FASTQ_FILE = "read1.fq"
READ2_FASTQ_FILE = "read2.fq"
FASTQ_FILES = (READ1_FASTQ_FILE, READ2_FASTQ_FILE)

# copy uncompress sequence step
COPY_UNCOMPRESS_JOB_MEM = 1000
COPY_UNCOMPRESS_JOB_WALLTIME = "10:00:00"

# fastqc
FASTQC_DIR_EXTENSION = "_fastqc"
FASTQC_DATA_FILE = "fastqc_data.txt"
FASTQC_REPORT_FILE = "fastqc_report.html"
FASTQC_JOB_MEM = 2000
FASTQC_JOB_WALLTIME = "10:00:00"

# abundant sequence mapping
ABUNDANT_SAM_FILES = ('abundant_hits_read1.sam', 'abundant_hits_read2.sam')
ABUNDANT_BAM_FILE = 'abundant_hits.bam'
SORTED_ABUNDANT_BAM_FILE = 'abundant_hits.srt.bam'
ABUNDANT_MAPPING_JOB_MEM = 2000
ABUNDANT_MAPPING_JOB_WALLTIME = "20:00:00"
ABUNDANT_FILTER_JOB_MEM = 1000
ABUNDANT_FILTER_JOB_WALLTIME = "10:00:00"
ABUNDANT_BAMSORT_JOB_MEM = 8192
ABUNDANT_BAMSORT_JOB_WALLTIME = "10:00:00"

# filtered fastq files
FILTERED_FASTQ_PREFIX = 'filtered_read'
FILTERED_FASTQ_FILES = tuple(("%s%d.fq" % (FILTERED_FASTQ_PREFIX,x)) for x in (1,2)) 

# contamination sequence mapping
XENO_SAM_FILES = ('xeno_hits_read1.sam', 'xeno_hits_read2.sam')
# filtered fastq files
XENO_BAM_FILE = 'xeno_hits.bam'
SORTED_XENO_BAM_FILE = 'xeno_hits.srt.bam'
XENO_MAPPING_JOB_MEM = 3750
XENO_MAPPING_JOB_WALLTIME = "20:00:00"
XENO_FILTER_JOB_MEM = 1000
XENO_FILTER_JOB_WALLTIME = "10:00:00"
XENO_BAMSORT_JOB_MEM = 8192
XENO_BAMSORT_JOB_WALLTIME = "10:00:00"

# fragment size distribution
FRAG_SIZE_DIST_FILE = "frag_size_dist.txt"
FRAG_SIZE_DIST_PLOT_FILE = "frag_size_dist_plot.pdf"
FRAG_SIZE_JOB_MEM = 3750
FRAG_SIZE_JOB_WALLTIME = "10:00:00"

# TODO: remove chimerascan from pipeline
# chimerascan gene fusion analysis
#CHIMERASCAN_DIR = 'chimerascan'
#CHIMERASCAN_RESULTS_FILE = 'chimeras.bedpe'
#CHIMERASCAN_MIN_SEGMENT_LENGTH = 25

# tophat alignment results
TOPHAT_DIR = 'tophat'
TOPHAT_BAM_FILE = os.path.join(TOPHAT_DIR, "accepted_hits.bam")
TOPHAT_JUNCTIONS_FILE = os.path.join(TOPHAT_DIR, "junctions.bed")
TOPHAT_JOB_MEM = 16000
TOPHAT_JOB_WALLTIME = "72:00:00"

# picard metrics files
LANE_ALIGNMENT_SUMMARY_METRICS = "lane.alignment_summary_metrics"
LANE_INSERT_SIZE_HISTOGRAM_PDF = "lane.insert_size_histogram.pdf"
LANE_INSERT_SIZE_METRICS = "lane.insert_size_metrics"
LANE_QUALITY_BY_CYCLE_METRICS = "lane.quality_by_cycle_metrics"
LANE_QUALITY_BY_CYCLE_PDF = "lane.quality_by_cycle.pdf"
LANE_QUALITY_DISTRIBUTION_METRICS = "lane.quality_distribution_metrics"
LANE_QUALITY_DISTRIBUTION_PDF = "lane.quality_distribution.pdf"
LANE_RNASEQ_METRICS = "lane.rnaseq_metrics"
LANE_RNASEQ_METRICS_PLOT_PDF = "lane.rnaseq_metrics_plot.pdf"
PICARD_METRICS_JOB_MEM = 8192
PICARD_METRICS_JOB_WALLTIME = "20:00:00"

# coverage bedgraph file
COVERAGE_BEDGRAPH_FILE = "coverage.bedgraph"
COVERAGE_BIGWIG_FILE = "coverage.bigwig"
COVERAGE_JOB_MEM = 3750
COVERAGE_JOB_WALLTIME = "20:00:00"
BIGWIG_JOB_MEM = 3750
BIGWIG_JOB_WALLTIME = "20:00:00"

# merged alignment results
MERGED_FRAG_SIZE_DIST_FILE = "merged_frag_size_dist.txt"
MERGED_BAM_FILE = "merged_alignments.bam"
MERGED_CLEAN_BAM_FILE = "merged_alignments_srdup.bam"
MERGE_FRAG_SIZE_JOB_MEM = 1000
MERGE_FRAG_SIZE_JOB_WALLTIME = "1:00:00"
MERGE_BAM_JOB_MEM = 1000
MERGE_BAM_JOB_WALLTIME = "1:00:00"
CLEAN_BAM_JOB_MEM = 3750
CLEAN_BAM_JOB_WALLTIME = "60:00:00"

# snp calling results
SAMTOOLS_VARIANT_BCF_FILE = "samtools.var.raw.bcf"
SAMTOOLS_VARIANT_VCF_FILE = "samtools.var.flt.vcf"
SAMTOOLS_VARIANT_JOB_MEM = 12000
SAMTOOLS_VARIANT_JOB_WALLTIME = "60:00:00"
VARSCAN_VARIANT_SNV_FILE = "varscan.snvs.flt.txt"
VARSCAN_VARIANT_IND_FILE = "varscan.indels.txt"
VARSCAN_VARIANT_JOB_MEM = 12000
VARSCAN_VARIANT_JOB_WALLTIME = "60:00:00"
BWA_ALIGNMENT_JOB_MEM=3750
BWA_ALIGNMENT_JOB_WALLTIME="24:00:00"
BWA_MAPPING_JOB_MEM=8192
BWA_MAPPING_JOB_WALLTIME='36:00:00'
SAM2BAM_JOB_MEM=12000
SAM2BAM_JOB_WALLTIME="24:00:00"
SORT_READS_JOB_MEM=12000
SORT_READS_JOB_WALLTIME="12:00:00"
MERGE_SAM_FILES_JOB_MEM=12000
MERGE_SAM_FILES_JOB_WALLTIME="24:00:00"
BAM_CLEANING_JOB_MEM=12000
BAM_CLEANING_JOB_WALLTIME="36:00:00"
SYMBOLIC_LINK_JOB_MEM=3750
SYMBOLIC_LINK_JOB_WALLTIME="1:00:00"
COSMIC_COVERAGE_JOB_MEM=3750
COSMIC_COVERAGE_JOB_WALLTIME="24:00:00"
CAPTURE_COVERAGE_JOB_MEM=3750
CAPTURE_COVERAGE_JOB_WALLTIME="24:00:00"
EXOME_CNV_JOB_MEM=3750
EXOME_CNV_JOB_WALLTIME="24:00:00"
TUMOR_COSMIC_VCF = "tumor_cosmic_positions.vcf"
COSMIC_QUAL_VCF = "cov_cosmic_positions.vcf"

# cufflinks output
CUFFLINKS_DIR = "cufflinks"
CUFFLINKS_TRANSCRIPTS_GTF_FILE = os.path.join(CUFFLINKS_DIR, "transcripts.gtf")
CUFFLINKS_GENES_FILE = os.path.join(CUFFLINKS_DIR, "genes.fpkm_tracking")
CUFFLINKS_ISOFORMS_FILE = os.path.join(CUFFLINKS_DIR, "isoforms.fpkm_tracking")
CUFFLINKS_JOB_MEM = 24000
CUFFLINKS_JOB_WALLTIME = "60:00:00"

# notify complete
NOTIFY_COMPLETE_JOB_MEM = 500
NOTIFY_COMPLETE_JOB_WALLTIME = "1:00:00"

# cleanup intermediate files
CLEANUP_INTERMEDIATE_FILES_JOB_MEM = 500
CLEANUP_INTERMEDIATE_FILES_JOB_WALLTIME = "1:00:00"

# cnv results
CNV_FILE = "exome.cnvs.txt"
LOH_FILE = "exome_loh.txt"
CNV_PLOT = "exome_cnv_plot"

# files for exome analysis
# alignment results
ALIGNED_READS_SAM = "aligned_reads.sam"
ALIGNED_READS_BAM = "aligned_reads.bam"
ALIGNED_READS_BAM_TMP = "aligned_reads_tmp.bam"
ALIGNED_BAM_SORTED = "aligned_reads_sorted.bam"
# Merged BAM FILES
MERGED_BAM_EFILE = "merged_aligned.bam"
MERGED_BAM_MDUP_EFILE = "merged_aligned_mdup.bam"
MERGED_CLEAN_BAM_EFILE = "merged_aligned_clean.bam"
# coverage
EXOME_COVERAGE = "coverage_exome.cov"
PROBE_COVERAGE = "coverage_probes.cov"
PROBE_COVERAGE_SUMMARY = "probes_coverage_summary.cov"

# job complete
JOB_COMPLETE_FILE = "job.done"
DNA_JOB_COMPLETE_FILE = "dnajobs.done"
RNA_JOB_COMPLETE_FILE = "rnajobs.done"

def get_tophat_library_type(strand_protocol):
    if strand_protocol == "dutp":
        return "fr-firststrand"
    else:
        return "fr-unstranded"

def get_picard_strand_specificity(strand_protocol):
    if strand_protocol == "dutp":
        return "SECOND_READ_TRANSCRIPTION_STRAND"
    else:
        return "NONE"

def attach_exome_library_to_results(library, root_dir):
    library.output_dir = os.path.join(root_dir, library.id)
    # merged alignment files
    library.merged_bam_efile = os.path.join(library.output_dir, MERGED_BAM_EFILE)   
    library.merged_bam_mdup_efile= os.path.join(library.output_dir, MERGED_BAM_MDUP_EFILE)
    library.merged_cleaned_bam_efile = os.path.join(library.output_dir, MERGED_CLEAN_BAM_EFILE)
    # lane results
    for lane in library.lanes:
        lane.output_dir = os.path.join(library.output_dir, lane.id)
        # FASTQ files
        lane.fastq_files = [lane.read1_file]
        lane.copied_fastq_files = [os.path.join(lane.output_dir, READ1_FASTQ_FILE)]
        if lane.read2_file is not None:
            lane.fastq_files.append(lane.read2_file)
            lane.copied_fastq_files.append(os.path.join(lane.output_dir, READ2_FASTQ_FILE))
        # FASTQC results
        lane.fastqc_data_files = []
        lane.fastqc_report_files = []
        for readnum in xrange(len(lane.copied_fastq_files)):
            fastqc_dir = os.path.join(lane.output_dir, "%s%s" % (os.path.basename(lane.fastq_files[readnum]), FASTQC_DIR_EXTENSION))
            lane.fastqc_data_files.append(os.path.join(fastqc_dir, FASTQC_DATA_FILE))
            lane.fastqc_report_files.append(os.path.join(fastqc_dir, FASTQC_REPORT_FILE))
        # EXOME alignment files
        lane.exome_sam_aln = os.path.join(lane.output_dir, ALIGNED_READS_SAM)
        lane.exome_bam_aln = os.path.join(lane.output_dir, ALIGNED_READS_BAM)
        lane.exome_bam_tmp = os.path.join(lane.output_dir, ALIGNED_READS_BAM_TMP)
        lane.exome_bam_sorted = os.path.join(lane.output_dir, ALIGNED_BAM_SORTED)

def attach_exome_sample_to_results(sample, root_dir):    
    sample.output_dir = os.path.join(root_dir, sample.id)
    sample.job_complete_file = os.path.join(sample.output_dir, JOB_COMPLETE_FILE)
    sample.xml_file = os.path.join(sample.output_dir, SAMPLE_XML_FILE)
    # output files produced by exome pipeline
    sample.merged_bam_efile = os.path.join(sample.output_dir, MERGED_BAM_EFILE)   
    sample.merged_bam_mdup_efile= os.path.join(sample.output_dir, MERGED_BAM_MDUP_EFILE)
    sample.merged_cleaned_bam_efile = os.path.join(sample.output_dir, MERGED_CLEAN_BAM_EFILE)
    sample.cosmic_qual_vcf_file = os.path.join(sample.output_dir, COSMIC_QUAL_VCF)
    sample.exome_coverage_file = os.path.join(sample.output_dir, EXOME_COVERAGE)
    sample.probe_coverage_file = os.path.join(sample.output_dir, PROBE_COVERAGE)
    sample.probe_summary_file = os.path.join(sample.output_dir, PROBE_COVERAGE_SUMMARY)
    sample.coverage_bedgraph_file = os.path.join(sample.output_dir, COVERAGE_BEDGRAPH_FILE)
    sample.coverage_bigwig_file = os.path.join(sample.output_dir, COVERAGE_BIGWIG_FILE)      
    # attach libraries to results
    for library in sample.libraries:
        attach_exome_library_to_results(library, sample.output_dir)


def attach_rnaseq_library_to_results(library, root_dir):
    library.output_dir = os.path.join(root_dir, library.id)
    # merged fragment size distribution
    library.merged_frag_size_dist_file = os.path.join(library.output_dir, MERGED_FRAG_SIZE_DIST_FILE)
    # merged BAM file
    library.merged_bam_file = os.path.join(library.output_dir, MERGED_BAM_FILE)   
    library.merged_cleaned_bam_file = os.path.join(library.output_dir, MERGED_CLEAN_BAM_FILE)   
    # rnaseq SNP calling files
    library.samtools_bcf_file = os.path.join(library.output_dir, SAMTOOLS_VARIANT_BCF_FILE)
    library.samtools_vcf_file = os.path.join(library.output_dir, SAMTOOLS_VARIANT_VCF_FILE)
    library.varscan_snv_file = os.path.join(library.output_dir, VARSCAN_VARIANT_SNV_FILE)
    library.varscan_indel_file = os.path.join(library.output_dir, VARSCAN_VARIANT_IND_FILE)    
    # Cufflinks output files
    library.cufflinks_dir = os.path.join(library.output_dir, CUFFLINKS_DIR)
    library.cufflinks_gtf_file = os.path.join(library.output_dir, CUFFLINKS_TRANSCRIPTS_GTF_FILE)
    library.cufflinks_genes_fpkm_file = os.path.join(library.output_dir, CUFFLINKS_GENES_FILE)
    library.cufflinks_isoforms_fpkm_file = os.path.join(library.output_dir, CUFFLINKS_ISOFORMS_FILE)  
    # lane results
    for lane in library.lanes:
        lane.output_dir = os.path.join(library.output_dir, lane.id)
        # FASTQ files
        lane.fastq_files = [lane.read1_file]
        lane.copied_fastq_files = [os.path.join(lane.output_dir, READ1_FASTQ_FILE)]
        if lane.read2_file is not None:
            lane.fastq_files.append(lane.read2_file)
            lane.copied_fastq_files.append(os.path.join(lane.output_dir, READ2_FASTQ_FILE))
        # FASTQC results
        lane.fastqc_data_files = []
        lane.fastqc_report_files = []
        for readnum in xrange(len(lane.fastq_files)):
            fastqc_dir = os.path.join(lane.output_dir, "%s%s" % (os.path.basename(lane.fastq_files[readnum]), FASTQC_DIR_EXTENSION))
            lane.fastqc_data_files.append(os.path.join(fastqc_dir, FASTQC_DATA_FILE))
            lane.fastqc_report_files.append(os.path.join(fastqc_dir, FASTQC_REPORT_FILE))
        # Abundant SAM files
        lane.abundant_sam_files = []
        for readnum in xrange(len(lane.copied_fastq_files)):
            lane.abundant_sam_files.append(os.path.join(lane.output_dir, ABUNDANT_SAM_FILES[readnum]))
        # Filtered abundant BAM and FASTQ
        lane.abundant_bam_file = os.path.join(lane.output_dir, ABUNDANT_BAM_FILE)
        lane.filtered_fastq_files = []
        for readnum in xrange(len(lane.copied_fastq_files)):
            lane.filtered_fastq_files.append(os.path.join(lane.output_dir, FILTERED_FASTQ_FILES[readnum]))
        # TODO: remove chimerascan
        # chimerascan results
        #lane.chimerascan_dir = os.path.join(lane.output_dir, CHIMERASCAN_DIR)
        #lane.chimerascan_results_file = os.path.join(lane.chimerascan_dir, CHIMERASCAN_RESULTS_FILE)
        # Sorted abundant reads bam file
        lane.sorted_abundant_bam_file = os.path.join(lane.output_dir, SORTED_ABUNDANT_BAM_FILE)
        # Contaminant foreign organism (xeno) SAM files
        lane.xeno_sam_files = []
        for readnum in xrange(len(lane.copied_fastq_files)):
            lane.xeno_sam_files.append(os.path.join(lane.output_dir, XENO_SAM_FILES[readnum]))
        # Contaminant foreign organism BAM files
        lane.xeno_bam_file = os.path.join(lane.output_dir, XENO_BAM_FILE)
        lane.sorted_xeno_bam_file = os.path.join(lane.output_dir, SORTED_XENO_BAM_FILE)
        # Fragment size distribution
        lane.frag_size_dist_file = os.path.join(lane.output_dir, FRAG_SIZE_DIST_FILE)
        lane.frag_size_dist_plot_file = os.path.join(lane.output_dir, FRAG_SIZE_DIST_PLOT_FILE)
        # Align reads with tophat
        lane.tophat_dir = os.path.join(lane.output_dir, TOPHAT_DIR)
        lane.tophat_bam_file = os.path.join(lane.output_dir, TOPHAT_BAM_FILE)
        lane.tophat_juncs_file = os.path.join(lane.output_dir, TOPHAT_JUNCTIONS_FILE)
        # Picard metrics
        lane.alignment_summary_metrics = os.path.join(lane.output_dir, LANE_ALIGNMENT_SUMMARY_METRICS)
        lane.insert_size_histogram_pdf = os.path.join(lane.output_dir, LANE_INSERT_SIZE_HISTOGRAM_PDF)
        lane.insert_size_metrics = os.path.join(lane.output_dir, LANE_INSERT_SIZE_METRICS)
        lane.quality_by_cycle_metrics = os.path.join(lane.output_dir, LANE_QUALITY_BY_CYCLE_METRICS)
        lane.quality_by_cycle_pdf = os.path.join(lane.output_dir, LANE_QUALITY_BY_CYCLE_PDF)
        lane.quality_distribution_metrics = os.path.join(lane.output_dir, LANE_QUALITY_DISTRIBUTION_METRICS)
        lane.quality_distribution_pdf = os.path.join(lane.output_dir, LANE_QUALITY_DISTRIBUTION_PDF)
        lane.rnaseq_metrics = os.path.join(lane.output_dir, LANE_RNASEQ_METRICS)
        lane.rnaseq_metrics_pdf = os.path.join(lane.output_dir, LANE_RNASEQ_METRICS_PLOT_PDF)
        # Coverage file
        lane.coverage_bedgraph_file = os.path.join(lane.output_dir, COVERAGE_BEDGRAPH_FILE)
        # Bigwig file
        lane.coverage_bigwig_file = os.path.join(lane.output_dir, COVERAGE_BIGWIG_FILE)

def attach_rnaseq_sample_to_results(sample, root_dir):    
    sample.output_dir = os.path.join(root_dir, sample.id)
    sample.job_complete_file = os.path.join(sample.output_dir, JOB_COMPLETE_FILE)
    sample.xml_file = os.path.join(sample.output_dir, SAMPLE_XML_FILE)
    # attach libraries to results
    for library in sample.libraries:
        attach_rnaseq_library_to_results(library, sample.output_dir)

def attach_sample_group_to_results(grp, root_dir):
    grp.output_dir = os.path.join(root_dir, grp.id)
    grp.xml_file = os.path.join(grp.output_dir, SAMPLE_GROUP_XML_FILE)
    grp.job_complete_file = os.path.join(grp.output_dir, JOB_COMPLETE_FILE)
    # output files produced by rna pipeline
    grp.rna_job_complete_file = os.path.join(grp.output_dir, RNA_JOB_COMPLETE_FILE)
    # output files produced by exome pipeline
    grp.dna_job_complete_file = os.path.join(grp.output_dir, DNA_JOB_COMPLETE_FILE)
    grp.samtools_bcf_file = os.path.join(grp.output_dir, SAMTOOLS_VARIANT_BCF_FILE)
    grp.samtools_vcf_file = os.path.join(grp.output_dir, SAMTOOLS_VARIANT_VCF_FILE)
    grp.varscan_snv_file = os.path.join(grp.output_dir, VARSCAN_VARIANT_SNV_FILE)
    grp.varscan_indel_file = os.path.join(grp.output_dir, VARSCAN_VARIANT_IND_FILE)        
    grp.tumor_cosmic_file = os.path.join(grp.output_dir, TUMOR_COSMIC_VCF)
    grp.exome_cnv_file = os.path.join(grp.output_dir, CNV_FILE)
    grp.exome_loh_file = os.path.join(grp.output_dir, LOH_FILE)
    grp.exome_cnv_plot = os.path.join(grp.output_dir, CNV_PLOT)
    # attach samples to results
    for sample_type, sample in grp.samples.iteritems():
        if sample is None:
            continue
        if ((sample_type == SAMPLE_TYPE_EXOME_TUMOR) or
            (sample_type == SAMPLE_TYPE_EXOME_NORMAL)):
            attach_exome_sample_to_results(sample, grp.output_dir)
        else:
            attach_rnaseq_sample_to_results(sample, grp.output_dir)

def attach_patient_to_results(patient, root_dir):
    patient.output_dir = os.path.join(root_dir, patient.id)
    patient.job_complete_file = os.path.join(patient.output_dir, JOB_COMPLETE_FILE)
    patient.xml_file = os.path.join(patient.output_dir, PATIENT_XML_FILE)
    for grp in patient.sample_groups.itervalues():
        attach_sample_group_to_results(grp, patient.output_dir)

def validate_patient_results(patient):
    is_valid = True
    for grp in patient.sample_groups.itervalues():
        for sample_type, sample in grp.samples.iteritems():
            if ((sample_type == SAMPLE_TYPE_RNASEQ) or
                (sample_type == SAMPLE_TYPE_CAPTURE_RNASEQ)):
                is_valid = is_valid and validate_rnaseq.validate_sample_results(sample)
    return is_valid

class AnalysisConfig(object):

    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        c = AnalysisConfig()        
        c.patients = []
        # read patients
        for patient_elem in root.findall("patient"):
            c.patients.append(Patient.from_xml(patient_elem))
        return c
    
    def to_xml(self, output_file):
        root = etree.Element("analysis")
        for patient in self.patients:
            patient.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        # write to file
        f = open(output_file, "w")
        print >>f, etree.tostring(root)
        f.close()            
    
    def is_valid(self):
        valid = True
        for patient in self.patients:
            valid = valid and patient.is_valid()
        return valid
    
    def attach_to_results(self, root_dir):
        for patient in self.patients:            
            attach_patient_to_results(patient, root_dir)
            
    def validate_results(self):
        for patient in self.patients:            
            validate_patient_results(patient)

class GenomeConfig(object):
    __fields__ = ("root_dir",
                  "abundant_bowtie_index",
                  "xeno_bowtie_index",
                  "genome_bowtie_index",
                  "genome_bowtie2_index",
                  "genome_fasta_file",
                  "genome_bwa_index",
                  "fragment_size_bowtie_index",
                  "gene_annotation_refflat",
                  "picard_ribosomal_intervals",
                  "chrom_sizes",
                  "cosmic_positions",
                  "capture_roche",
                  "capture_agilent",
                  "exome_bed_file")
    
    @staticmethod
    def from_xml_elem(elem):
        g = GenomeConfig()
        g.species = elem.get("name")
        for attrname in GenomeConfig.__fields__:
            setattr(g, attrname, elem.findtext(attrname))
        return g

    def to_xml(self, root):
        root.set("name", self.species)
        for attrname in GenomeConfig.__fields__:
            elem = etree.SubElement(root, attrname)
            elem.text = str(getattr(self, attrname))            

    def get_path(self, attr_name):
        return os.path.join(self.root_dir, getattr(self, attr_name))
    
    def is_valid(self, references_dir=""):
        valid = True
        abs_root_dir = os.path.join(references_dir, self.root_dir)
        if not os.path.exists(abs_root_dir):
            logging.error("genome root directory %s not found" % (self.root_dir))
            valid = False            
        if not os.path.exists(os.path.join(abs_root_dir, self.abundant_bowtie_index + ".1.bt2")):
            logging.error("Abundant bowtie2 index %s not found" % (self.abundant_bowtie_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.abundant_bowtie_index + ".fa")):
            logging.error("Abundant sequences fasta file %s not found" % (self.abundant_bowtie_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.xeno_bowtie_index + ".1.bt2")):
            logging.error("Foreign contaminants bowtie2 index %s not found" % (self.xeno_bowtie_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.genome_bowtie_index + ".1.ebwt")):
            logging.error("Genome bowtie index %s not found" % (self.genome_bowtie_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.genome_bowtie2_index + ".1.bt2")):
            logging.error("Genome bowtie2 index %s not found" % (self.genome_bowtie2_index))
            valid = False
        if not os.path.exists(os.path.join(abs_root_dir, self.fragment_size_bowtie_index + ".1.ebwt")):
            logging.error("Fragment size bowtie index %s not found" % (self.fragment_size_bowtie_index))
            valid = False
        if not os.path.join(abs_root_dir, self.genome_bwa_index + ".bwt"):
            logging.error("Genome BWA index %s not found" % (self.genome_bwa_index))
            valid = False
        for attrname in ("genome_fasta_file",
                         "gene_annotation_refflat",
                         "picard_ribosomal_intervals",
                         "chrom_sizes"):
            if not os.path.exists(os.path.join(abs_root_dir, getattr(self, attrname))):
                logging.error("Genome file %s not found" % (getattr(self, attrname)))
                valid = False
        return valid

class ServerConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        c = ServerConfig()
        c.name = elem.get("name")
        c.address = elem.get("address", None)
        c.ssh_port = int(elem.get("ssh_port", "22"))
        c.modules_init_script = elem.findtext("modules_init_script")
        c.output_dir = elem.findtext("output_dir")
        c.tmp_dir = elem.findtext("tmp_dir")
        c.seq_repo_mirror_dir = elem.findtext("seq_repo_mirror_dir", None)
        c.references_dir = elem.findtext("references_dir")
        c.node_mem = float(elem.findtext("node_mem"))
        c.node_processors = int(elem.findtext("node_processors"))
        has_pbs = elem.findtext("pbs").lower()
        if has_pbs == "yes":
            c.pbs = True
            c.max_user_jobs = int(elem.findtext("max_user_jobs"))
        else:
            c.pbs = False
            c.max_user_jobs = None
        c.pbs_script_lines = []
        for line_elem in elem.findall("script_line"):
            c.pbs_script_lines.append(line_elem.text)
        return c
    
    def to_xml(self, root):
        root.set("name", self.name)
        root.set("address", str(self.address))
        root.set("ssh_port", str(self.ssh_port))
        for attrname in ("modules_init_script",
                         "output_dir",
                         "tmp_dir",
                         "references_dir",
                         "seq_repo_mirror_dir",
                         "node_mem",
                         "node_processors"):
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        if self.pbs:
            elem = etree.SubElement(root, "pbs")
            elem.text = "yes"
            elem = etree.SubElement(root, "max_user_jobs")
            elem.text = str(self.max_user_jobs)
        else:
            elem = etree.SubElement(root, "pbs")
            elem.text = "no"
        for line in self.pbs_script_lines:
            elem = etree.SubElement(root, "script_line")
            elem.text = line

    def is_valid(self):
        valid = True
        # check directories
        if os.path.isdir(self.tmp_dir):
            logging.debug("Checking for 'tmp' directory '%s'... found" % (self.tmp_dir))
        else:
            logging.error("'tmp' directory '%s' not found" % (self.tmp_dir))
            valid = False
        if os.path.isdir(self.output_dir):
            logging.debug("Checking for 'output' directory %s... found" % (self.output_dir))
        else:
            logging.error("'output' directory %s not found" % (self.output_dir))
            valid = False            
        if os.path.isdir(self.references_dir):
            logging.debug("Checking for 'references' directory %s... found" % (self.references_dir))
        else:
            logging.error("'references' directory %s not found" % (self.references_dir))
            valid = False
        if os.path.exists(self.modules_init_script):
            logging.debug("Checking for modules init script %s... found" % (self.modules_init_script))
        else:
            logging.error("modules init script %s not found" % (self.modules_init_script))
            valid = False
        return valid    

# TODO: remove chimerascan?
#class ChimerascanConfig(object):
#    @staticmethod
#    def from_xml_elem(elem):
#        c = ChimerascanConfig()
#        for attrname in ("index", "trim5", "trim3", "frag_size_percentile"):
#            setattr(c, attrname, elem.findtext(attrname))
#        c.args = []
#        for arg_elem in elem.findall("arg"):
#            c.args.append(arg_elem.text)
#        return c
#    
#    def to_xml(self, root):
#        for attrname in ("index", "trim5", "trim3", "frag_size_percentile"):
#            attrval = getattr(self,attrname)
#            if attrval is not None:
#                elem = etree.SubElement(root, attrname)
#                elem.text = str(attrval)
#        for arg in self.args:
#            elem = etree.SubElement(root, "arg")
#            elem.text = arg            
#
#    def is_valid(self, species_dir):
#        valid = True
#        species_sub = lambda arg: arg.replace("${SPECIES}", species_dir)
#        newindex = species_sub(self.index)
#        if not os.path.exists(newindex):
#            logging.error("File not found: %s" % (newindex))
#            valid = False
#        return valid

class BWAConfig(object):
    __fields__ = ("bwa_cores", "nmismatch", "bwa_qual_trim", "mapping_qual")

    @staticmethod
    def from_xml_elem(elem):
        c = BWAConfig()
        for attrname in BWAConfig.__fields__:
            setattr(c, attrname, elem.findtext(attrname))
        c.args = []
        for arg_elem in elem.findall("arg"):
            c.args.append(arg_elem.text)
        return c

    def to_xml(self, root):
        for attrname in BWAConfig.__fields__:
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for arg in self.args:
            elem = etree.SubElement(root, "arg")
            elem.text = arg
                      
    def is_valid(self):
        # TODO: write this function
        return True


class VarscanConfig(object):
    __fields__ = ("min_coverage", "min_reads_alt", "min_avgbase_quality", 
                  "min_var_freq", "min_pvalue")
    @staticmethod
    def from_xml_elem(elem):
        c = VarscanConfig()
        for attrname in VarscanConfig.__fields__:
            setattr(c, attrname, elem.findtext(attrname))
        c.args = []
        for arg_elem in elem.findall("arg"):
            c.args.append(arg_elem.text)
        return c

    def to_xml(self, root):
        for attrname in VarscanConfig.__fields__:
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for arg in self.args:
            elem = etree.SubElement(root, "arg")
            elem.text = arg  

    def is_valid(self):
        # TODO: write this function
        return True

        
class PipelineConfig(object):
    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)  
        root = tree.getroot()
        c = PipelineConfig()
        # modules
        modules_elem = root.find("modules")
        c.modules = []
        for elem in modules_elem.findall("module"):
            c.modules.append(elem.text)
        # binaries (hard-coded for now)
        c.fastqc_bin = "fastqc"
        c.samtools_bin = "samtools"
        c.r_bin = "R"
        c.rscript_bin = "Rscript"
        c.bowtie_bin = "bowtie"
        c.bowtie2_bin = "bowtie2"
        c.tophat_bin = "tophat"
        c.cufflinks_bin = "cufflinks"
        c.bedtools_dir = ""
        c.ucsc_dir = ""
        c.picard_dir = ""
        c.varscan_dir = ""
        c.bwa_bin = "bwa"
        # check environment
        if "PICARDPATH" in os.environ:
            c.picard_dir = os.environ["PICARDPATH"]
        if "VARSCANPATH" in os.environ:
            c.varscan_dir = os.environ["VARSCANPATH"]
        # default fragment size parameters
        c.fragment_size_mean_default = int(root.findtext("fragment_size_mean_default"))
        c.fragment_size_stdev_default = int(root.findtext("fragment_size_stdev_default"))
        c.adaptor_length_default = int(root.findtext("adaptor_length_default"))
        c.min_fragment_size = int(root.findtext("min_fragment_size"))
        c.max_fragment_size = int(root.findtext("max_fragment_size"))
        # tophat parameters
        c.tophat_args = []
        elem = root.find("tophat")
        for arg_elem in elem.findall("arg"):
            c.tophat_args.append(arg_elem.text)
        # cufflinks parameters
        c.cufflinks_args = []
        elem = root.find("cufflinks")
        for arg_elem in elem.findall("arg"):
            c.cufflinks_args.append(arg_elem.text)
        # TODO: chimerascan has now been deprecated, remove this
        # chimerascan parameters
        # c.chimerascan_config = ChimerascanConfig.from_xml_elem(root.find("chimerascan"))
        # server setup
        c.servers = {}
        for elem in root.findall("server"):
            server = ServerConfig.from_xml_elem(elem)
            c.servers[server.name] = server
        # bwa parameters 
        c.bwa_config = BWAConfig.from_xml_elem(root.find("bwa"))
        # varscan parameters
        c.vscan_config = VarscanConfig.from_xml_elem(root.find("varscan"))
        # genome config
        c.species = {}
        for elem in root.findall("species"):
            g = GenomeConfig.from_xml_elem(elem)
            c.species[g.species] = g
        return c

    def to_xml(self, output_file):
        root = etree.Element("pipeline")
        # modules
        modules_elem = etree.SubElement(root, "modules")
        for m in self.modules:
            elem = etree.SubElement(modules_elem, "module")
            elem.text = m
        # fragment size parameters
        for attrname in ("fragment_size_mean_default",
                         "fragment_size_stdev_default",
                         "adaptor_length_default",
                         "min_fragment_size",
                         "max_fragment_size"):
            elem = etree.SubElement(root, attrname)
            elem.text = str(getattr(self, attrname))
        # tophat parameters
        tophat_elem = etree.SubElement(root, "tophat")
        for arg in self.tophat_args:
            elem = etree.SubElement(tophat_elem, "arg")
            elem.text = arg
        # cufflinks parameters
        cufflinks_elem = etree.SubElement(root, "cufflinks")
        for arg in self.cufflinks_args:
            elem = etree.SubElement(cufflinks_elem, "arg")
            elem.text = arg
        # TODO: chimerascan has now been deprecated, remove this
        # chimerascan parameters
        #chimerascan_elem = etree.SubElement(root, "chimerascan")
        #self.chimerascan_config.to_xml(chimerascan_elem)
        # bwa parameters
        bwa_elem = etree.SubElement(root, "bwa")
        self.bwa_config.to_xml(bwa_elem)
        # varscan parameters
        vscan_elem = etree.SubElement(root, "varscan")
        self.vscan_config.to_xml(vscan_elem)
        # servers
        for server in self.servers.itervalues():            
            elem = etree.SubElement(root, "server")
            server.to_xml(elem)
        # genomes
        for genome in self.species.itervalues():
            elem = etree.SubElement(root, "species")
            genome.to_xml(elem)
        # output files
        f = open(output_file, "w")
        # indent for pretty printing
        indent_xml(root)
        print >>f, etree.tostring(root)
        f.close() 

    def is_valid(self, server_name):
        """
        ensure configuration is valid
        """
        valid = True
        # check server
        if server_name not in self.servers:
            logging.error("Server %s not found" % (server_name))
            return False
        server = self.servers[server_name]
        if not server.is_valid():
            logging.error("Server %s missing required paths" % (server_name))
            valid = False
        # check genomes
        for species_name,genome in self.species.iteritems():
            if not genome.is_valid(server.references_dir):
                logging.error("Genome %s missing required files" % (species_name))
                valid = False
        # TODO: chimerascan has now been deprecated, remove this
        # check chimerascan config
        #species_dir = os.path.join(server.references_dir, genome.root_dir)
        #if not self.chimerascan_config.is_valid(species_dir):
        #    logging.error("Chimerascan missing required files")
        #    valid = False
        # check bwa config
        if not self.bwa_config.is_valid():
            logging.error("Invalid BWA configuration")
            valid = False
        # check varscan config
        if not self.vscan_config.is_valid():
            logging.error("Invalid Varscan configuration")
            valid = False
        #
        # Check software installation
        #
        # check fastqc
        msg = 'fastqc'
        if check_executable(self.fastqc_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # check samtools
        msg = 'samtools'
        if check_executable(self.samtools_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # java
        msg = 'java'
        if check_executable("java"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # picard
        jarfile = os.path.join(self.picard_dir, "CollectMultipleMetrics.jar")
        if os.path.exists(jarfile):
            logging.debug("Checking for picard tools... found")
        else:
            logging.error("Picard jarfile '%s' not found" % (jarfile))
            valid = False
        # varscan
        jarfile = os.path.join(self.varscan_dir, "varscan")
        if os.path.exists(jarfile):
            logging.debug("Checking for varscan... found")
        else:
            logging.error("varscan jarfile '%s' not found" % (jarfile))
            valid = False
        # check R
        msg = 'R'
        if check_executable(self.r_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'Rscript'
        if check_executable(self.rscript_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False                
        msg = 'bowtie'
        if check_executable(self.bowtie_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'bowtie2'
        if check_executable(self.bowtie2_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'tophat'
        if check_executable(self.tophat_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        msg = 'cufflinks'
        if check_executable(self.cufflinks_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # check BEDTools
        msg = 'BEDTools'
        if check_executable(os.path.join(self.bedtools_dir, "intersectBed")):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False            
        # check UCSC binaries
        msg = 'UCSC bedGraphToBigWig'
        if check_executable(os.path.join(self.ucsc_dir, "bedGraphToBigWig")):
            logging.debug("Checking for '%s' binaries... found" % msg)
        else:
            logging.error("'%s' binaries not found or not executable" % msg)
            valid = False
        msg = 'UCSC blat'
        if check_executable(os.path.join(self.ucsc_dir, "blat")):
            logging.debug("Checking for '%s' binaries... found" % msg)
        else:
            logging.error("'%s' binaries not found or not executable" % msg)
            valid = False
        msg = 'UCSC faToTwoBit'
        if check_executable(os.path.join(self.ucsc_dir, "faToTwoBit")):
            logging.debug("Checking for '%s' binaries... found" % msg)
        else:
            logging.error("'%s' binaries not found or not executable" % msg)
            valid = False
        # check for bx python library
        try:
            import bx.intervals.intersection
            logging.debug("Checking for 'bx python' library... found")
        except ImportError, e:
            logging.error("Package 'bx python' not found")
            valid = False
        # check for pysam library
        try:
            import pysam
            logging.debug("Checking for 'pysam' library... found")
        except ImportError, e:
            logging.error("Package 'pysam' not found")
            valid = False
        # TODO: chimerascan deprecated so remove this
        # chimerascan binary
        #msg = 'chimerascan'
        #if check_executable("chimerascan_run.py"):
        #    logging.debug("Checking for '%s' binary... found" % msg)
        #else:
        #    logging.error("'%s' binary not found or not executable" % msg)
        #    valid = False            
        # check for chimerascan libraries
        #try:
        #    import chimerascan
        #    logging.debug("Checking for 'chimerascan' library... found")
        #except ImportError, e:
        #    logging.error("Package 'chimerascan' not found")
        #    valid = False
        # check bwa
        msg = 'BWA'
        if check_executable(self.bwa_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
        # TODO: check for varscan executable
        return valid
