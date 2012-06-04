'''
Created on Aug 16, 2011

@author: mkiyer
'''
import argparse
import logging
import sys
import os
import xml.etree.cElementTree as etree

from oncoseq.lib import config
from oncoseq.lib.seqdb import Sample
from oncoseq.lib.config import attach_rnaseq_sample_to_results
from oncoseq.rnaseq.lib.validate import validate_sample_results

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_xml_file")
    parser.add_argument("output_dir")
    parser.add_argument("sample_id")
    args = parser.parse_args()
    # read samples
    tree = etree.parse(args.sample_xml_file)        
    root = tree.getroot()
    result = config.JOB_SUCCESS
    for elem in root.findall("sample"):
        sample = Sample.from_xml(elem)
        attach_rnaseq_sample_to_results(sample, args.output_dir)   
        # cleanup sample
        if sample.id != args.sample_id:
            continue
        if not validate_sample_results(sample):
            logging.error("Sample %s not valid - skipping cleanup step" % (sample.id))
            result = config.JOB_ERROR
            continue
        for library in sample.libraries:
            for lane in library.lanes:
                # remove fastq input files                 
                for f in lane.copied_fastq_files:
                    if os.path.exists(f):
                        os.remove(f)
                for f in lane.filtered_fastq_files:
                    if os.path.exists(f):
                        os.remove(f)
                # remove abundant SAM/BAM files (leave only the sorted BAM file)
                for f in lane.abundant_sam_files:
                    if os.path.exists(f):
                        os.remove(f)
                if os.path.exists(lane.abundant_bam_file):
                    os.remove(lane.abundant_bam_file)
                # remote foreign contam SAM/BAM files (leave only the sorted BAM file)
                for f in lane.xeno_sam_files:
                    if os.path.exists(f):
                        os.remove(f)
                if os.path.exists(lane.xeno_bam_file):
                    os.remove(lane.xeno_bam_file)                
                # remove coverage bedgraph file
                if os.path.exists(lane.coverage_bedgraph_file):
                    os.remove(lane.coverage_bedgraph_file)
    return result
    
if __name__ == '__main__':
    sys.exit(main())