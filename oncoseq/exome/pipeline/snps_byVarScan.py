'''
Created on Feb 15, 2012

@author: oabalbin
'''
import logging
import argparse
import subprocess
import sys
import os

from vs2vcf import convert_varscan_to_vcf
from oncoseq.lib import config
#from oncoseq.lib.vs2vcf import read_file

def call_snps(ref_fa, tumor_bam_file, 
              snps_vcf_file, indels_vcf_file,varscan_dir):
    # call snps and store in BCF file
    '''
    '''
    
    # TODO: Improve on using standard names, and do not use .replace() to change the extension of different files
    snps_vcf_file_tmp = snps_vcf_file.replace('.flt.txt','.raw.txt')
    min_mapq=str(20)
    min_baseq=str(20)
    mpileup_t=tumor_bam_file.replace('.bam','.mpileup.bcf')
    mpileup_f=open(mpileup_t,'w')
    #mpileup_f=open(bcf_file.replace('.bcf','.mpileup.bcf'),'w')
    args = ["samtools", "mpileup", "-q",min_mapq,"-Q",min_baseq,
            "-f",ref_fa, tumor_bam_file] #+['>', mpileup_f]
    mpileup_p = subprocess.call(args, stdout=mpileup_f)
    mpileup_f.close()
    print args
    if mpileup_p !=0:
        logging.error("samtools mpileup tumor failed")
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
            return config.JOB_ERROR  
        
    # Call SNVs using varscan
    dpval,tumor_fraction= 0.10,0.05 # default pvalue to force pvalue calculation
    #f=open(snps_vcf_file,'w')
    #f.close()
    args=["java","-jar",os.path.join(varscan_dir,"varscan"),"mpileup2snp", mpileup_t,
          "--variants",str(1),
          "--p-value",str(dpval),
          "--output-vcf", str(0)]    
    print args
    f = open(snps_vcf_file_tmp, "w")    
    retcode = subprocess.call(args,stdout=f)
    f.close()
    
    if retcode != 0:
        logging.error("varscan SNV calling failed")
        if os.path.exists(snps_vcf_file):
            os.remove(snps_vcf_file)
            return config.JOB_ERROR
    
    args=["java","-jar",os.path.join(varscan_dir,"varscan"),"mpileup2indel", mpileup_t,
          "--variants",str(1),
          "--p-value",str(dpval),
          "--output-vcf", str(0)]

    print args
    f = open(indels_vcf_file, "w")    
    retcode = subprocess.call(args,stdout=f)
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
    parser.add_argument("--varscan-dir", dest="varscan_dir", default="")
    parser.add_argument("ref_fa")
    parser.add_argument("tumor_bam_file")
    parser.add_argument("bcf_file")
    parser.add_argument("vcf_file")
    args = parser.parse_args()
    return call_snps(args.ref_fa, args.tumor_bam_file, 
                     args.bcf_file, args.vcf_file,varscan_dir=args.varscan_dir)


if __name__ == '__main__': 
    sys.exit(main())