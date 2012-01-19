'''
Created on Aug 14, 2011

@author: mkiyer
'''
import logging
import argparse
import subprocess
import sys
import os

from oncoseq.lib import config

def call_snps(ref_fa, bam_file, bcf_file, vcf_file):
    # call snps and store in BCF file
    args = ["samtools", "mpileup", "-uf",
            ref_fa, bam_file]
    mpileup_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    args = ["bcftools", "view", "-bvcg", "-"]
    f = open(bcf_file, "w")
    retcode = subprocess.call(args, stdin=mpileup_p.stdout, stdout=f)
    f.close()
    if retcode != 0:
        logging.error("samtools snp caller failed")
        mpileup_p.kill()
        if os.path.exists(bcf_file):
            os.remove(bcf_file)
            return config.JOB_ERROR            
    # convert to VCF file
    args = ["bcftools", "view", bcf_file]
    bcfview_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    args = ["vcfutils.pl", "varFilter"]
    f = open(vcf_file, "w")
    retcode = subprocess.call(args, stdin=bcfview_p.stdout, stdout=f)
    f.close()
    if retcode != 0:
        logging.error("samtools VCF generation failed")
        bcfview_p.kill()
        if os.path.exists(vcf_file):
            os.remove(vcf_file)
            return config.JOB_ERROR            
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("ref_fa")
    parser.add_argument("bam_file")
    parser.add_argument("bcf_file")
    parser.add_argument("vcf_file")
    args = parser.parse_args()
    return call_snps(args.ref_fa, args.bam_file, args.bcf_file, args.vcf_file)

if __name__ == '__main__': 
    sys.exit(main())