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

def call_snps(ref_fa, normal_bam_file, tumor_bam_file, 
              snps_vcf_file, indels_vcf_file,varscan_dir):
    # call snps and store in BCF file
    '''
    '''
    
    # TODO: Improve on using standard names, and do not use .replace() to change the extension of different files
    snps_vcf_file_tmp = snps_vcf_file.replace('.flt.txt','.raw.txt')
    min_mapq=str(20)
    min_baseq=str(20)
    mpileup_t=tumor_bam_file.replace('.bam','.mpileup.bcf')
    mpileup_n=normal_bam_file.replace('.bam','.mpileup.bcf')
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
    
    mpileup_f=open(mpileup_n,'w')
    #mpileup_f=open(bcf_file.replace('.bcf','.mpileup.bcf'),'w')
    args = ["samtools", "mpileup", "-q",min_mapq,"-Q",min_baseq,
            "-f",ref_fa, normal_bam_file] #+['>', mpileup_f]
    mpileup_p = subprocess.call(args, stdout=mpileup_f)
    mpileup_f.close()
    print args
    if mpileup_p !=0:
        logging.error("samtools mpileup normal failed")
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
            return config.JOB_ERROR  
    
    
    # Call SNVs using varscan
    dpval,tumor_fraction= 0.10,0.05 # default pvalue to force pvalue calculation
    #f=open(snps_vcf_file,'w')
    #f.close()
    args=["java","-jar",os.path.join(varscan_dir,"varscan"),"somatic", mpileup_n,mpileup_t,  
          "--variants",str(1),
          "--p-value",str(dpval),
          "--output-snp",snps_vcf_file_tmp,
          "--output-indel", indels_vcf_file]#,'>',snps_vcf_file] #"--output-vcf",str(1) for output in vcf format
    args = ",".join(args).replace(',',' ')
    print args
    retcode = subprocess.call(args, shell=True)
    
    if retcode != 0:
        logging.error("varscan SNV calling failed")
        if os.path.exists(snps_vcf_file):
            os.remove(snps_vcf_file)
            return config.JOB_ERROR
    
    # Filter somatic SNVs using varscan
    args=["java","-jar",os.path.join(varscan_dir,"varscan"),"somaticFilter",
          snps_vcf_file_tmp,
          "--indel-file",indels_vcf_file,
          "--output-file",snps_vcf_file,
          "--p-value",str(dpval),
          "--min-var-freq",str(tumor_fraction)]
    args = ",".join(args).replace(',',' ')
    print args
    retcode = subprocess.call(args,shell=True)
    if retcode != 0:
        logging.error("varscan Somatic Filter calling failed")
        if os.path.exists(indels_vcf_file) and os.path.exists(snps_vcf_file):
            os.remove(indels_vcf_file)
            os.remove(snps_vcf_file)
            return config.JOB_ERROR
    else:
        convert_varscan_to_vcf(snps_vcf_file,snps_vcf_file.replace('.txt','.vcf'))
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
            os.remove(mpileup_n)
    
    
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--varscan-dir", dest="varscan_dir", default="")
    parser.add_argument("ref_fa")
    parser.add_argument("normal_bam_file")
    parser.add_argument("tumor_bam_file")
    parser.add_argument("bcf_file")
    parser.add_argument("vcf_file")
    args = parser.parse_args()
    return call_snps(args.ref_fa, args.normal_bam_file, args.tumor_bam_file, 
                     args.bcf_file, args.vcf_file,varscan_dir=args.varscan_dir)


if __name__ == '__main__': 
    sys.exit(main())