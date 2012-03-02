'''
Created on Feb 1, 2012

@author: oabalbin
'''
import logging
import argparse
import subprocess
import sys
import os

from oncoseq.lib import config
from oncoseq.lib.vs2vcf import read_file

def call_snps(ref_fa, bam_file, snps_vcf_file, indels_vcf_file):
    # call snps and store in BCF file
    min_mapq=str(20)
    min_baseq=str(20)
    mpileup_t=bam_file.replace('.bam','.mpileup.bcf')
    mpileup_f=open(mpileup_t,'w')
    #mpileup_f=open(bcf_file.replace('.bcf','.mpileup.bcf'),'w')
    args = ["samtools", "mpileup", "-q",min_mapq,"-Q",min_baseq,
            "-f",ref_fa, bam_file] #+['>', mpileup_f]
    
    mpileup_p = subprocess.call(args, stdout=mpileup_f)
    mpileup_f.close()
    
    if mpileup_p !=0:
        logging.error("samtools mpileup failed")
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
            return config.JOB_ERROR  
    
    # Call SNVs using varscan
    dpval= 0.10 # default pvalue to force pvalue calculation
    #f=open(snps_vcf_file,'w')
    #f.close()
    args=["java","-jar",os.path.join("$VARSCANPATH","varscan"),"mpileup2snp", mpileup_t,   
          "--variants",str(1),
          "--p-value",str(dpval)]#,'>',snps_vcf_file] #"--output-vcf",str(1) for output in vcf format
    args = ",".join(args).replace(',',' ')
    print args
    f=open(snps_vcf_file,'w')
    retcode = subprocess.call(args, stdout=f,shell=True)
    f.close()
    if retcode != 0:
        logging.error("varscan SNV calling failed")
        if os.path.exists(snps_vcf_file):
            os.remove(snps_vcf_file)
            return config.JOB_ERROR
    else:
        read_file(snps_vcf_file,snps_vcf_file+'.vcf')

    #f=open(indels_vcf_file,'w')
    #f.close()        
    args=["java","-jar",os.path.join("$VARSCANPATH","varscan"),"mpileup2indel",mpileup_t,
          "--variants",str(1),
          "--p-value",str(dpval),
          "--output-vcf",str(1)]#,'>',indels_vcf_file] #"--output-vcf",str(1) for output in vcf format
    args = ",".join(args).replace(',',' ')
    print args
    f=open(indels_vcf_file,'w')
    retcode = subprocess.call(args, stdout=f,shell=True)
    f.close()
    if retcode != 0:
        logging.error("varscan Indel calling failed")
        if os.path.exists(indels_vcf_file):
            os.remove(indels_vcf_file)
            return config.JOB_ERROR
    else:
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
    
    
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