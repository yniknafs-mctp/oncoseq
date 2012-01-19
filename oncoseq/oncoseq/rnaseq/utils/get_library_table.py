'''
Created on Nov 28, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging

from oncoseq.lib.seqdb import SeqDB
from oncoseq.lib.base import file_exists_and_nz_size
from oncoseq.lib.config import attach_sample_to_results, \
    validate_library_results, check_sam_file, get_tophat_library_type, \
    check_frag_size_dist_file, FRAGMENT_LAYOUT_PAIRED
from oncoseq.lib.fragment_size_distribution import FragmentSizeDistribution
        
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("seqdb_xls_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    seqdb = SeqDB.from_xls(args.seqdb_xls_file)
    # print header
    print '\t'.join(["cohort", "patient", "sample", "library", "lanes",
                     "library_type", "frag_len_mean", "frag_len_std_dev",
                     "has_pe_lanes", "valid", "cufflinks_gtf_file", 
                     "bam_file"])    
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
            # get library type
            #
            library_type = get_tophat_library_type(library.strand_protocol)
            #
            # find whether any of the lanes are paired-end
            #
            has_pe_lanes = any((lane.fragment_layout == FRAGMENT_LAYOUT_PAIRED) 
                               for lane in library.lanes)
            #
            # get fragment size distribution parameters
            # 
            frag_len_mean = "NA"
            frag_len_std_dev = "NA"
            if check_frag_size_dist_file(library.merged_frag_size_dist_file):
                # get fragment size parameters
                frag_size_dist = FragmentSizeDistribution.from_file(open(library.merged_frag_size_dist_file))
                frag_len_mean = str(int(frag_size_dist.isize_at_percentile(50.0)))
                frag_len_std_dev = str(int(round(frag_size_dist.std(),0)))
            #
            # get merged bam file
            #
            if check_sam_file(library.merged_bam_file, isbam=True):
                bam_file = library.merged_bam_file
            else:
                logging.warning("Sample %s Library %s missing BAM file" % (sample.id, library.id))
                bam_file = "NA"
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
            # cohort, patient, sample, library, lane, is_valid, aligned_reads, read_length, junctions_file
            fields = [sample.cohort, sample.patient_id, sample.id, library.id,
                      ",".join([lane.id for lane in library.lanes]), 
                      library_type, frag_len_mean, frag_len_std_dev, 
                      int(has_pe_lanes),
                      int(is_valid),
                      cufflinks_gtf_file,
                      bam_file]
            print '\t'.join(map(str, fields))
                
if __name__ == '__main__': 
    sys.exit(main())
