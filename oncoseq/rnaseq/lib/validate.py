'''
Created on May 21, 2012

@author: mkiyer
'''
import logging

from oncoseq.lib.fragment_size_distribution import FragmentSizeDistribution
from oncoseq.lib.seqdb import FRAGMENT_LAYOUT_PAIRED
from oncoseq.lib.base import file_exists_and_nz_size
from oncoseq.lib.sam import check_sam_file

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
    # TODO: remove chimerascan
    # check chimerascan results (only run chimerascan for paired-end reads)    
    #if ((len(lane.filtered_fastq_files) > 1) and 
    #    (not file_exists_and_nz_size(lane.chimerascan_results_file))):
    #    logging.error("Lane %s missing/corrupt chimerascan results file" % (lane.id))
    #    is_valid = False
    # check sorted foreign sequence bam file
    #if not check_sam_file(lane.sorted_xeno_bam_file, isbam=True):
    # TODO: 09-13-2012 tmp.
    if not check_sam_file(lane.xeno_bam_file, isbam=True):
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
    # check snp files from varscan
    if not file_exists_and_nz_size(library.varscan_snv_file):
        logging.error("Library %s missing varscan snv file" % (library.id))
        is_valid = False
    return is_valid

def validate_sample_results(sample):
    is_valid = True
    for library in sample.libraries:
        is_valid = is_valid and validate_library_results(library)
    return is_valid