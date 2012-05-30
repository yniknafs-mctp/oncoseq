'''
Created on Aug 3, 2011

@author: mkiyer
'''
import os
import logging
import xml.etree.cElementTree as etree
from collections import defaultdict
from base import check_executable, indent_xml, file_exists_and_nz_size
from seqdb import Patient, Sample, Library, Lane, FRAGMENT_LAYOUT_PAIRED
from fragment_size_distribution import FragmentSizeDistribution

# job return codes
JOB_SUCCESS = 0
JOB_ERROR = 1

# Protocols/analysis Folders
EXOME = "exome"
EXOME_DIR="exome"
RNASEQ = "rnaseq"
RNASEQ_DIR="rnaseq"
RNASEQ_CAPTURE="capture_rnaseq"
RNASEQ_CAPTURE_DIR="capture_rnaseq"

# Protocols supported in oncoseq 
VALID_PROTOCOLS={EXOME:EXOME_DIR,RNASEQ:RNASEQ_DIR,
                 RNASEQ_CAPTURE:RNASEQ_CAPTURE_DIR}


# remote job constants
REMOTE_ANALYSIS_XML_FILE = "analysis.xml"
REMOTE_CONFIG_XML_FILE = "config.xml"
REMOTE_CODE_TARGZ_FILE = "code.tar.gz"

# job file constants
READ1_FASTQ_FILE = "read1.fq"
READ2_FASTQ_FILE = "read2.fq"
FASTQ_FILES = (READ1_FASTQ_FILE, READ2_FASTQ_FILE)

# fastqc
FASTQC_READ1_DIR = "read1.fq_fastqc"
FASTQC_READ1_DATA_FILE = os.path.join(FASTQC_READ1_DIR, 'fastqc_data.txt')
FASTQC_READ1_REPORT_FILE = os.path.join(FASTQC_READ1_DIR, 'fastqc_report.html')
FASTQC_READ2_DIR = "read2.fq_fastqc"
FASTQC_READ2_DATA_FILE = os.path.join(FASTQC_READ2_DIR, 'fastqc_data.txt')
FASTQC_READ2_REPORT_FILE = os.path.join(FASTQC_READ2_DIR, 'fastqc_report.html')
FASTQC_DATA_FILES = [FASTQC_READ1_DATA_FILE, FASTQC_READ2_DATA_FILE]
FASTQC_REPORT_FILES = [FASTQC_READ1_REPORT_FILE, FASTQC_READ2_REPORT_FILE]

# abundant sequence mapping
ABUNDANT_SAM_FILES = ('abundant_hits_read1.sam', 'abundant_hits_read2.sam')
# filtered fastq files
ABUNDANT_BAM_FILE = 'abundant_hits.bam'
SORTED_ABUNDANT_BAM_FILE = 'abundant_hits.srt.bam'
FILTERED_FASTQ_PREFIX = 'filtered_read'
FILTERED_FASTQ_FILES = tuple(("%s%d.fq" % (FILTERED_FASTQ_PREFIX,x)) for x in (1,2)) 

# contamination sequence mapping
XENO_SAM_FILES = ('xeno_hits_read1.sam', 'xeno_hits_read2.sam')
# filtered fastq files
XENO_BAM_FILE = 'xeno_hits.bam'
SORTED_XENO_BAM_FILE = 'xeno_hits.srt.bam'

# fragment size distribution
FRAG_SIZE_DIST_FILE = "frag_size_dist.txt"
FRAG_SIZE_DIST_PLOT_FILE = "frag_size_dist_plot.pdf"

# chimerascan gene fusion analysis
CHIMERASCAN_DIR = 'chimerascan'
CHIMERASCAN_RESULTS_FILE = 'chimeras.bedpe'
CHIMERASCAN_MIN_SEGMENT_LENGTH = 25

# tophat alignment results
TOPHAT_DIR = 'tophat'
TOPHAT_BAM_FILE = os.path.join(TOPHAT_DIR, "accepted_hits.bam")
TOPHAT_JUNCTIONS_FILE = os.path.join(TOPHAT_DIR, "junctions.bed")

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

# coverage bedgraph file
COVERAGE_BEDGRAPH_FILE = "coverage.bedgraph"
COVERAGE_BIGWIG_FILE = "coverage.bigwig"

# merged alignment results
MERGED_BAM_FILE = "merged_alignments.bam"
MERGED_CLEAN_BAM_FILE = "merged_alignments_srdup.bam"
MERGED_FRAG_SIZE_DIST_FILE = "merged_frag_size_dist.txt"

# snp calling results
SAMTOOLS_VARIANT_BCF_FILE = "samtools.var.raw.bcf"
SAMTOOLS_VARIANT_VCF_FILE = "samtools.var.flt.vcf"
VARSCAN_VARIANT_SNV_FILE = "varscan.snvs.flt.txt"
VARSCAN_VARIANT_IND_FILE = "varscan.indels.txt"
TUMOR_COSMIC_VCF = "tumor_cosmic_positions.vcf"
COSMIC_QUAL_VCF = "cov_cosmic_positions.vcf"

# CNV results
CNV_FILE = "exome.cnvs.txt"
LOH_FILE = "exome_loh.txt"
CNV_PLOT = "exome_cnv_plot"

# cufflinks output
CUFFLINKS_DIR = "cufflinks"
CUFFLINKS_TRANSCRIPTS_GTF_FILE = os.path.join(CUFFLINKS_DIR, "transcripts.gtf")
CUFFLINKS_GENES_FILE = os.path.join(CUFFLINKS_DIR, "genes.fpkm_tracking")
CUFFLINKS_ISOFORMS_FILE = os.path.join(CUFFLINKS_DIR, "isoforms.fpkm_tracking")


# Files for EXOME analysis
# alignment results
ALIGNED_READS_SAM="aligned_reads.sam"
ALIGNED_READS_BAM="aligned_reads.bam"
ALIGNED_READS_BAM_TMP="aligned_reads_tmp.bam"
ALIGNED_BAM_SORTED="aligned_reads_sorted.bam"
# Merged BAM FILES
MERGED_BAM_EFILE="merged_aligned.bam"
MERGED_BAM_MDUP_EFILE="merged_aligned_mdup.bam"
MERGED_CLEAN_BAM_EFILE="merged_aligned_clean.bam"

# coverage
EXOME_COVERAGE="coverage_exome.cov"
PROBE_COVERAGE="coverage_probes.cov"
PROBE_COVERAGE_SUMMARY="probes_coverage_summary.cov"
# job complete
JOB_COMPLETE_FILE = "job.done"

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


def attach_patient_to_results(patient, samples, root_dir):
    patient.output_dir = os.path.join(root_dir, patient.id)
    patient.exome_dir=os.path.join(patient.output_dir,EXOME_DIR)
    patient.rnaseq_dir=os.path.join(patient.output_dir,RNASEQ_DIR)
    patient.capture_rnaseq_dir=os.path.join(patient.output_dir,RNASEQ_CAPTURE_DIR)
    
    # SNV output files
    patient.job_complete_file = os.path.join(patient.output_dir, JOB_COMPLETE_FILE)
    patient.samtools_bcf_file = os.path.join(patient.exome_dir, SAMTOOLS_VARIANT_BCF_FILE)
    patient.samtools_vcf_file = os.path.join(patient.exome_dir, SAMTOOLS_VARIANT_VCF_FILE)
    patient.varscan_snv_file = os.path.join(patient.exome_dir, VARSCAN_VARIANT_SNV_FILE)
    patient.varscan_indel_file = os.path.join(patient.exome_dir, VARSCAN_VARIANT_IND_FILE)        
    patient.tumor_cosmic_file = os.path.join(patient.exome_dir,TUMOR_COSMIC_VCF)
    patient.exome_cnv_file = os.path.join(patient.exome_dir,CNV_FILE)
    patient.exome_loh_file = os.path.join(patient.exome_dir,LOH_FILE)
    patient.exome_cnv_plot = os.path.join(patient.exome_dir,CNV_PLOT)

    # Include GATK files
    
    # attach samples to results
    attach_sample_to_results(samples, patient.output_dir)
    

def attach_sample_to_results(samples, root_dir):        
    
    for sample in samples:
        
        protocol_path = os.path.join(root_dir,VALID_PROTOCOLS[sample.protocol])
        sample.output_dir = os.path.join(protocol_path, sample.id)
        sample.job_complete_file = os.path.join(sample.output_dir, JOB_COMPLETE_FILE)

        # Files for EXOME output
        sample.merged_bam_efile = os.path.join(sample.output_dir, MERGED_BAM_EFILE)   
        sample.merged_bam_mdup_efile= os.path.join(sample.output_dir, MERGED_BAM_MDUP_EFILE)
        sample.merged_cleaned_bam_efile = os.path.join(sample.output_dir, MERGED_CLEAN_BAM_EFILE)
        sample.cosmic_qual_vcf_file = os.path.join(sample.output_dir, COSMIC_QUAL_VCF)
        sample.exome_coverage_file = os.path.join(sample.output_dir,EXOME_COVERAGE)
        sample.probe_coverage_file = os.path.join(sample.output_dir,PROBE_COVERAGE)
        sample.probe_summary_file = os.path.join(sample.output_dir,PROBE_COVERAGE_SUMMARY)
        sample.coverage_bedgraph_file = os.path.join(sample.output_dir,COVERAGE_BEDGRAPH_FILE)
        sample.coverage_bigwig_file = os.path.join(sample.output_dir,COVERAGE_BIGWIG_FILE)
        
        # attach libraries to results
        attach_library_to_results(sample, sample.output_dir)
        

def attach_library_to_results(sample, root_dir):
        
    for library in sample.libraries:
        library.output_dir = os.path.join(sample.output_dir, library.id)
        # merged fragment size distribution
        library.merged_frag_size_dist_file = os.path.join(library.output_dir, MERGED_FRAG_SIZE_DIST_FILE)
        # merged BAM file
        library.merged_bam_file = os.path.join(library.output_dir, MERGED_BAM_FILE)   
        library.merged_cleaned_bam_file = os.path.join(library.output_dir, MERGED_CLEAN_BAM_FILE)
        
        # Files for EXOME output
        library.merged_bam_efile = os.path.join(library.output_dir, MERGED_BAM_EFILE)   
        library.merged_bam_mdup_efile= os.path.join(library.output_dir, MERGED_BAM_MDUP_EFILE)
        library.merged_cleaned_bam_efile = os.path.join(library.output_dir, MERGED_CLEAN_BAM_EFILE)

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
            for readnum in xrange(len(lane.copied_fastq_files)):
                lane.fastqc_data_files.append(os.path.join(lane.output_dir, FASTQC_DATA_FILES[readnum]))
                lane.fastqc_report_files.append(os.path.join(lane.output_dir, FASTQC_REPORT_FILES[readnum]))
            # Abundant SAM files
            lane.abundant_sam_files = []
            for readnum in xrange(len(lane.copied_fastq_files)):
                lane.abundant_sam_files.append(os.path.join(lane.output_dir, ABUNDANT_SAM_FILES[readnum]))
            # Filtered abundant BAM and FASTQ
            lane.abundant_bam_file = os.path.join(lane.output_dir, ABUNDANT_BAM_FILE)
            lane.filtered_fastq_files = []
            for readnum in xrange(len(lane.copied_fastq_files)):
                lane.filtered_fastq_files.append(os.path.join(lane.output_dir, FILTERED_FASTQ_FILES[readnum]))
            # chimerascan results
            lane.chimerascan_dir = os.path.join(lane.output_dir, CHIMERASCAN_DIR)
            lane.chimerascan_results_file = os.path.join(lane.chimerascan_dir, CHIMERASCAN_RESULTS_FILE)
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
            
            # EXOME alignment files
            lane.exome_sam_aln = os.path.join(lane.output_dir, ALIGNED_READS_SAM)
            lane.exome_bam_aln = os.path.join(lane.output_dir, ALIGNED_READS_BAM)
            lane.exome_bam_tmp = os.path.join(lane.output_dir, ALIGNED_READS_BAM_TMP)
            lane.exome_bam_sorted = os.path.join(lane.output_dir, ALIGNED_BAM_SORTED)

def check_tophat_juncs_file(filename):
    if not file_exists_and_nz_size(filename):
        return False
    lines = 0
    for line in open(filename):
        if not line or line.startswith("track"):
            continue
        lines += 1
        if lines > 100:
            return True
    return lines > 0

def check_sam_file(filename, isbam=False):
    is_valid = True
    if not file_exists_and_nz_size(filename):
        is_valid = False
    else:
        import pysam
        try:
            fmt = "rb" if isbam else "r"
            samfh = pysam.Samfile(filename, fmt)
            samfh.close()   
        except:
            is_valid = False
    return is_valid

def check_frag_size_dist_file(filename):
    is_valid=True
    if not file_exists_and_nz_size(filename):
        is_valid = False
    else:
        try:
            frag_size_dist = FragmentSizeDistribution.from_file(open(filename))
        except:
            is_valid = False    
    return is_valid

def validate_lane_results(lane):
    is_valid = True
    # check FASTQC results
    for f in lane.fastqc_data_files:
        if not file_exists_and_nz_size(f):
            logging.error("Lane %s missing fastqc data file %s" % (lane.id, f))
            is_valid = False
    for f in lane.fastqc_report_files:
        if not file_exists_and_nz_size(f):
            logging.error("Lane %s missing fastqc report file %s" % (lane.id, f))
            is_valid = False
    # check sorted abundant reads bam file
    if not check_sam_file(lane.sorted_abundant_bam_file, isbam=True):
        logging.error("Lane %s missing/corrupt abundant reads BAM file" % (lane.id))
        is_valid = False
    # check chimerascan results (only run chimerascan for paired-end reads)    
    if ((len(lane.filtered_fastq_files) > 1) and 
        (not file_exists_and_nz_size(lane.chimerascan_results_file))):
        logging.error("Lane %s missing/corrupt chimerascan results file" % (lane.id))
        is_valid = False
    # check sorted foreign sequence bam file
    if not check_sam_file(lane.sorted_xeno_bam_file, isbam=True):
        logging.error("Lane %s missing/corrupt foreign sequences BAM file" % (lane.id))
        is_valid = False
    # check fragment size distribution
    if not check_frag_size_dist_file(lane.frag_size_dist_file):
        logging.error("Lane %s missing/corrupt fragment size distribution" % (lane.id))
        is_valid = False
    # check tophat junctions file
    if not file_exists_and_nz_size(lane.tophat_juncs_file):
        logging.error("Lane %s missing/empty tophat junctions file" % (lane.id))
        is_valid = False
    # check tophat bam file
    if not check_sam_file(lane.tophat_bam_file, isbam=True):
        logging.error("Lane %s missing/corrupt tophat BAM file" % (lane.id))
        is_valid = False
    # check picard summary metrics
    if not file_exists_and_nz_size(lane.alignment_summary_metrics):
        logging.error("Lane %s missing picard alignment summary metrics" % (lane.id))
        is_valid = False
    if ((lane.fragment_layout == FRAGMENT_LAYOUT_PAIRED) and 
        (not file_exists_and_nz_size(lane.insert_size_metrics))):
        logging.error("Lane %s missing picard insert size metrics" % (lane.id))
        is_valid = False
    if not file_exists_and_nz_size(lane.quality_by_cycle_metrics):
        logging.error("Lane %s missing picard quality by cycle metrics" % (lane.id))
        is_valid = False
    if not file_exists_and_nz_size(lane.quality_distribution_metrics):
        logging.error("Lane %s missing picard quality distribution metrics" % (lane.id))
        is_valid = False
    if not file_exists_and_nz_size(lane.rnaseq_metrics):
        logging.error("Lane %s missing picard rnaseq metrics" % (lane.id))
        is_valid = False
    # check coverage files
    if not file_exists_and_nz_size(lane.coverage_bigwig_file):
        logging.error("Lane %s missing coverage bigwig file" % (lane.id))
        is_valid = False
    return is_valid

def validate_library_results(library):
    is_valid = True
    # validate child lanes
    for lane in library.lanes:
        is_valid = is_valid and validate_lane_results(lane)
    # check merged fragment size distribution file
    if not check_frag_size_dist_file(library.merged_frag_size_dist_file):
        logging.error("Library %s missing/corrupt fragment size distribution" % (library.id))
        is_valid = False
    # check merged BAM file
    if not check_sam_file(library.merged_bam_file, isbam=True):
        logging.error("Library %s missing/corrupt tophat BAM file" % (library.id))
        is_valid = False
    # check cufflinks files
    if not file_exists_and_nz_size(library.cufflinks_gtf_file):
        logging.error("Library %s missing cufflinks gtf file" % (library.id))
        is_valid = False
    # check snp files from samtools
    if not file_exists_and_nz_size(library.samtools_vcf_file):
        logging.error("Library %s missing samtools variant vcf file" % (library.id))
        is_valid = False
    # check snp files from varascan
    if not file_exists_and_nz_size(library.varscan_snv_file):
        logging.error("Library %s missing varscan snv file" % (library.id))
        is_valid = False
    return is_valid

def validate_sample_results(sample):
    is_valid = True
    for library in sample.libraries:
        is_valid = is_valid and validate_library_results(library)
    return is_valid

class AnalysisConfig(object):
    @staticmethod
    def from_xml(xmlfile):
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        c = AnalysisConfig()
        c.patients = defaultdict(list)
        c.samples = []
        c.species = set() 
        
        # read patients
        # TODO extend to multiples patients by analysis file. Usually one patient by analysis file
        for patient_elem in root.findall("patient"):
            patient = Patient.from_xml(patient_elem)
            c.species.add(patient.species)             
            #c.patients.append(patient)
            
            # read samples
            for sample_elem in patient_elem.findall("sample"):
                sample = Sample.from_xml(sample_elem)
                c.samples.append(sample)
                c.patients[patient].append(sample)
                
                for lib_elem in sample_elem.findall("library"):
                    lib = Library.from_xml(lib_elem)
                    lib.sample = sample
                    sample.libraries.append(lib)
                    for lane_elem in lib_elem.findall("lane"):
                        lane = Lane.from_xml(lane_elem)
                        lane.library = lib
                        lib.lanes.append(lane)
        
        if len(c.species) > 1:
            logging.error("Analysis include different species" % (patient.id))
        else:
            c.species=list(c.species)[0]
            
        return c
    
    def to_xml(self, output_file):
        root = etree.Element("analysis")
        # add patient species, gender age
        '''
        elem = etree.SubElement(root, "id")
        elem.text = self.patient.id
        elem = etree.SubElement(root, "species")
        elem.text = self.species
        elem = etree.SubElement(root, "gender")
        elem.text = self.gender
        elem = etree.SubElement(root, "age")
        elem.text = self.age
        '''
        # add samples
        for patient in self.patients:
            patient_elem = patient.to_xml(root)
            for sample in self.samples:
                sample_elem = sample.to_xml(patient_elem)
                for library in sample.libraries:
                    lib_elem = library.to_xml(sample_elem)
                    for lane in library.lanes:
                        lane.to_xml(lib_elem)
        f = open(output_file, "w")
        # indent for pretty printing
        indent_xml(root)
        print >>f, etree.tostring(root)
        f.close()            
    
    def is_valid(self):
        valid = True
        for sample in self.samples:
            valid = sample.is_valid() 
            for library in sample.libraries:        
                valid = library.is_valid() 
                for lane in library.lanes:
                    valid = lane.is_valid() 
        return valid
    
    def attach_to_results(self, root_dir):
        
        # Modified from the original folder structure
        for patient, samples in self.patients.iteritems():
            attach_patient_to_results(patient, samples, root_dir)
    

class GenomeConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        g = GenomeConfig()
        g.species = elem.get("name")
        for attrname in ("root_dir",
                         "abundant_bowtie_index",
                         "xeno_bowtie_index",
                         "genome_bowtie_index",
                         "genome_fasta_file",
                         "genome_bwa_index",
                         "fragment_size_bowtie_index",
                         "gene_annotation_refflat",
                         "picard_ribosomal_intervals",
                         "chrom_sizes",
                         "cosmic_positions",
                         "capture_roche",
                         "capture_agilent",
                         "exome_bed_file"):
            setattr(g, attrname, elem.findtext(attrname))
        return g

    def to_xml(self, root):
        root.set("name", self.species)
        for attrname in ("root_dir",
                         "abundant_bowtie_index",
                         "xeno_bowtie_index",
                         "genome_bowtie_index",
                         "genome_fasta_file",
                         "genome_bwa_index",
                         "fragment_size_bowtie_index",
                         "gene_annotation_refflat",
                         "picard_ribosomal_intervals",
                         "chrom_sizes",
                         "cosmic_positions",
                         "capture_roche",
                         "capture_agilent",
                         "exome_bed_file"):
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

class ChimerascanConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        c = ChimerascanConfig()
        for attrname in ("index", "trim5", "trim3", "frag_size_percentile"):
            setattr(c, attrname, elem.findtext(attrname))
        c.args = []
        for arg_elem in elem.findall("arg"):
            c.args.append(arg_elem.text)
        return c
    
    def to_xml(self, root):
        for attrname in ("index", "trim5", "trim3", "frag_size_percentile"):
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for arg in self.args:
            elem = etree.SubElement(root, "arg")
            elem.text = arg            

    def is_valid(self, species_dir):
        valid = True
        species_sub = lambda arg: arg.replace("${SPECIES}", species_dir)
        newindex = species_sub(self.index)
        if not os.path.exists(newindex):
            logging.error("File not found: %s" % (newindex))
            valid = False
        return valid   

class BWAConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        c = BWAConfig()
        for attrname in ("bwa_cores", "nmismatch", "bwa_qual_trim", "mapping_qual"):
            setattr(c, attrname, elem.findtext(attrname))
        c.args = []
        for arg_elem in elem.findall("arg"):
            c.args.append(arg_elem.text)
        return c

    def to_xml(self, root):
        for attrname in ("bwa_cores", "nmismatch", "bwa_qual_trim", "mapping_qual"):
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for arg in self.args:
            elem = etree.SubElement(root, "arg")
            elem.text = arg  
                      
    # TODO: IS valid function
    #def is_valid(self, species_dir):


class VarScanConfig(object):
    @staticmethod
    def from_xml_elem(elem):
        c = VarScanConfig()
        for attrname in ("min_coverage", "min_reads_alt", "min_avgbase_quality", "min_var_freq","min_pvalue"):
            setattr(c, attrname, elem.findtext(attrname))
        c.args = []
        for arg_elem in elem.findall("arg"):
            c.args.append(arg_elem.text)
        return c

    def to_xml(self, root):
        for attrname in ("min_coverage", "min_reads_alt", "min_avgbase_quality", "min_var_freq","min_pvalue"):
            attrval = getattr(self,attrname)
            if attrval is not None:
                elem = etree.SubElement(root, attrname)
                elem.text = str(attrval)
        for arg in self.args:
            elem = etree.SubElement(root, "arg")
            elem.text = arg  
                      
    # TODO: IS valid function
    #def is_valid(self, species_dir):

        
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
        c.bwa_bin="bwa"
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
        # chimerascan parameters
        c.chimerascan_config = ChimerascanConfig.from_xml_elem(root.find("chimerascan"))
        # server setup
        c.servers = {}
        for elem in root.findall("server"):
            server = ServerConfig.from_xml_elem(elem)
            c.servers[server.name] = server
        
        # bwa parameters 
        c.bwa_config = BWAConfig.from_xml_elem(root.find("bwa"))
        # varscan parameters
        c.vscan_config = VarScanConfig.from_xml_elem(root.find("varscan"))
        
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
        # chimerascan parameters
        chimerascan_elem = etree.SubElement(root, "chimerascan")
        self.chimerascan_config.to_xml(chimerascan_elem)
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

    def is_valid(self, server_name, species_name):
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
        # check genome
        print species_name
        if species_name not in self.species:
            logging.error("Genome %s not found" % (species_name))
            return False
        genome = self.species[species_name]
        if not genome.is_valid(server.references_dir):
            logging.error("Genome %s missing required files" % (species_name))
            valid = False
        # check chimerascan config
        species_dir = os.path.join(server.references_dir, genome.root_dir)
        if not self.chimerascan_config.is_valid(species_dir):
            logging.error("Chimerascan missing required files")
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
        # chimerascan binary
        msg = 'chimerascan'
        if check_executable("chimerascan_run.py"):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False            
        # check for chimerascan libraries
        try:
            import chimerascan
            logging.debug("Checking for 'chimerascan' library... found")
        except ImportError, e:
            logging.error("Package 'chimerascan' not found")
            valid = False
        # check bwa
        msg = 'BWA'
        if check_executable(self.bwa_bin):
            logging.debug("Checking for '%s' binary... found" % msg)
        else:
            logging.error("'%s' binary not found or not executable" % msg)
            valid = False
    
            
        return valid
