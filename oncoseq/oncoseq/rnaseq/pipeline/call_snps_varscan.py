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
from oncoseq.lib.vs2vcf import convert_varscan_to_vcf

def call_snps(ref_fa, bam_file, 
              varscan_snv_file, 
              varscan_indel_file,
              varscan_dir):
    #
    # Use samtools mpileup to create a BCF file of snps
    #
    min_mapq = "20"
    min_baseq = "20"
    mpileup_t = bam_file.replace('.bam','.mpileup.bcf')
    mpileup_f = open(mpileup_t,'w')
    #mpileup_f=open(bcf_file.replace('.bcf','.mpileup.bcf'),'w')
    args = ["samtools", "mpileup", "-q", min_mapq,"-Q", min_baseq,
            "-f", ref_fa, bam_file] #+['>', mpileup_f]    
    logging.debug("samtools mpileup args: %s" % (args))
    retcode = subprocess.call(args, stdout=mpileup_f)
    mpileup_f.close()
    if retcode !=0:
        logging.error("samtools mpileup failed")
        if os.path.exists(mpileup_t):
            os.remove(mpileup_t)
            return config.JOB_ERROR  
    #
    # Call SNVs using varscan
    #
    # default pvalue to force pvalue calculation
    dpval = 0.10 
    args = ["java", "-jar", os.path.join(varscan_dir, "varscan"), 
            "mpileup2snp", mpileup_t, "--variants", "1", 
            "--p-value", str(dpval)]
    cmd = " ".join(map(str,args))
    logging.debug("VarScan args: %s" % (cmd))
    f = open(varscan_snv_file, 'w')
    retcode = subprocess.call(cmd, stdout=f, shell=True)
    f.close()
    if retcode != 0:
        logging.error("VarScan SNV calling failed")
        if os.path.exists(varscan_snv_file):
            os.remove(varscan_snv_file)
        return config.JOB_ERROR
    #
    # Convert the VarScan output file to VCF format
    # 
    convert_varscan_to_vcf(varscan_snv_file, varscan_snv_file + '.vcf')
    #
    # Call Indels using varscan
    #
    args = ["java","-jar",os.path.join("$VARSCANPATH","varscan"),"mpileup2indel",mpileup_t,
            "--variants", "1", 
            "--p-value", str(dpval)]
    #,'>',indels_vcf_file] #"--output-vcf",str(1) for output in vcf format
    cmd = " ".join(map(str,args))
    logging.debug("VarScan indel args: %s" % (cmd))
    f = open(varscan_indel_file,'w')
    retcode = subprocess.call(cmd, stdout=f,shell=True)
    f.close()
    if retcode != 0:
        logging.error("VarScan Indel calling failed")
        if os.path.exists(varscan_indel_file):
            os.remove(varscan_indel_file)
        return config.JOB_ERROR
    # remove intermediate files
    if os.path.exists(mpileup_t):
        os.remove(mpileup_t)
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--varscan-dir", dest="varscan_dir", default="")
    parser.add_argument("ref_fa")
    parser.add_argument("bam_file")
    parser.add_argument("snv_file")
    parser.add_argument("indel_file")
    args = parser.parse_args()
    return call_snps(args.ref_fa, args.bam_file, 
                     args.snv_file, args.indel_file,
                     varscan_dir=args.varscan_dir)

if __name__ == '__main__': 
    sys.exit(main())