'''
Created on Aug 16, 2011

@author: mkiyer
'''
import argparse
import logging
import sys
import os

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, validate_sample_results

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--rm-fastq", action="store_true", dest="rm_fastq", default=False)
    parser.add_argument("analysis_file")
    parser.add_argument("output_dir")
    parser.add_argument("sample_id")
    args = parser.parse_args()
    analysis = AnalysisConfig.from_xml(args.analysis_file)
    analysis.attach_to_results(args.output_dir)
    result = config.JOB_SUCCESS
    for sample in analysis.samples:
        if sample.id != args.sample_id:
            continue
        if not validate_sample_results(sample):
            logging.error("Sample %s not valid -- skipping cleanup step" % (sample.id))
            result = config.JOB_ERROR
            continue
        for library in sample.libraries:
            for lane in library.lanes:
                # remove fastq input files
                if args.rm_fastq:
                    for f in lane.fastq_files:
                        if os.path.exists(f):
                            os.remove(f)                    
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