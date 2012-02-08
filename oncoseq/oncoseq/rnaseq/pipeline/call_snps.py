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
    '''
    Varfilter default options:
    Options: -Q INT    minimum RMS mapping quality for SNPs [10]
         -d INT    minimum read depth [2]
         -D INT    maximum read depth [10000000]
         -a INT    minimum number of alternate bases [2]
         -w INT    SNP within INT bp around a gap to be filtered [3]
         -W INT    window size for filtering adjacent gaps [10]
         -1 FLOAT  min P-value for strand bias (given PV4) [0.0001]
         -2 FLOAT  min P-value for baseQ bias [1e-100]
         -3 FLOAT  min P-value for mapQ bias [0]
         -4 FLOAT  min P-value for end distance bias [0.0001]
                 -e FLOAT  min P-value for HWE (plus F<0) [0.0001]
         -p        print filtered variants
    '''

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