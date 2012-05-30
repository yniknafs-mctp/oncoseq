'''
Created on Nov 2, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging

from oncoseq.lib.seqdb import SeqDB
from oncoseq.lib.base import file_exists_and_nz_size
from oncoseq.lib.config import attach_sample_to_results, check_tophat_juncs_file, validate_lane_results, validate_library_results
from oncoseq.lib import picard
from oncoseq.lib import fastqc
from oncoseq.lib.seq import detect_read_length

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("seqdb_xls_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    seqdb = SeqDB.from_xls(args.seqdb_xls_file)
    # print header
    print '\t'.join(["cohort", "patient", "sample", "library", "lane", 
                     "valid", "aligned_reads", "read_length", 
                     "tophat_juncs_file", "cufflinks_gtf_file"])    
    for sample in seqdb.samples.itervalues():
        attach_sample_to_results(sample, os.path.abspath(args.output_dir))       
        for library in sample.libraries:
            #
            # check validity of library
            # 
            if validate_library_results(library):
                logging.info("Sample %s Library %s VALID" % (sample.id, library.id))
                is_valid = True
            else:
                logging.warning("Sample %s Library %s not valid" % (sample.id, library.id))
                is_valid = False
            #
            # get cufflinks GTF file 
            #
            if file_exists_and_nz_size(library.cufflinks_gtf_file):
                cufflinks_gtf_file = library.cufflinks_gtf_file
            else:
                logging.warning("Sample %s Library %s missing cufflinks GTF file" % (sample.id, library.id))
                cufflinks_gtf_file = "NA"
                is_valid = False            
            for lane in library.lanes:
                #
                # check validity of lane
                # 
                if lane.qc == "FAIL":
                    logging.info("Sample %s Lane %s failed QC" % (sample.id, lane.id))
                    continue
                # check validity of lane
                if validate_lane_results(lane):
                    logging.info("Sample %s Lane %s VALID" % (sample.id, lane.id))
                    is_valid = True
                else:
                    logging.warning("Sample %s Lane %s not valid" % (sample.id, lane.id))
                    is_valid = False
                #
                # get read length
                #
                if file_exists_and_nz_size(lane.fastqc_data_files[0]):
                    read_length = fastqc.get_read_length(lane.fastqc_data_files[0])
                else:
                    logging.error("Sample %s Lane %s missing fastqc data file" % (sample.id, lane.id))
                    if file_exists_and_nz_size(lane.fastq_files[0]):          
                        read_length = detect_read_length(lane.fastq_files[0])                    
                    else:
                        read_length = "NA"
                #
                # get number of aligned reads
                #
                if file_exists_and_nz_size(lane.alignment_summary_metrics):
                    num_aligned_reads = picard.get_num_aligned_reads(lane.alignment_summary_metrics)
                else:
                    logging.error("Sample %s Lane %s missing alignment summary metrics" % (sample.id, lane.id))
                    num_aligned_reads = "NA"
                #
                # get tophat junctions file
                #
                if check_tophat_juncs_file(lane.tophat_juncs_file):
                    tophat_juncs_file = lane.tophat_juncs_file
                else:
                    logging.error("Sample %s Lane %s missing tophat junctions file" % (sample.id, lane.id))
                    tophat_juncs_file = "NA"
                # cohort, patient, sample, library, lane, is_valid, aligned_reads, read_length, junctions_file
                fields = [sample.cohort, sample.patient_id, sample.id, 
                          library.id, lane.id, int(is_valid),
                          num_aligned_reads, 
                          read_length, 
                          tophat_juncs_file,
                          cufflinks_gtf_file]
                print '\t'.join(map(str, fields))
                
if __name__ == '__main__': 
    sys.exit(main())
