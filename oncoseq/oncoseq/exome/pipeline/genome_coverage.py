'''
Created on Oct 3, 2012

@author: oabalbin
'''


import logging
import argparse
import subprocess
import sys
import os
from oncoseq.lib import config

def genome_coverage(chrom_sizes,merged_cleaned_bam_efile,genomeCoverageFile):
    #It assumes the module bedtools is loaded in the node that is going to execute
    args = ["genomeCoverageBed",
            "-bg", "-split", "-g", chrom_sizes,
            "-ibam", merged_cleaned_bam_efile]
    
    f = open(genomeCoverageFile, "w")
    retcode = subprocess.call(args,stdout=f)
    f.close()
    
    if retcode != 0:
        logging.error("genomeCoverage computation failed")
        if os.path.exists(genomeCoverageFile):
            os.remove(genomeCoverageFile)
            return config.JOB_ERROR
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    #parser.add_argument("bedtools_dir")
    parser.add_argument("chrom_sizes")
    parser.add_argument("merged_cleaned_bam_efile")
    parser.add_argument("genomeCoverageFile")
    
    args = parser.parse_args()
    return genome_coverage(args.chrom_sizes,
                                args.merged_cleaned_bam_efile,
                                args.genomeCoverageFile)
        
    
if __name__ == '__main__':     
    sys.exit(main())