'''
Created on Apr 17, 2012

@author: oabalbin
'''
import logging
import argparse
import subprocess
import sys
import os

from oncoseq.lib import config
import oncoseq.lib
import oncoseq.exome.pipeline
import oncoseq.rnaseq.pipeline
_lib_dir=oncoseq.lib.__path__[0]
_pipeline_dir = oncoseq.exome.pipeline.__path__[0] 

def run_exomeCNV(samtools_bin, Rscrip_dir, benign_cov, tumor_cov, 
                     loh_file, cnv_file,cnv_plot,r_backup_file,
                     plot_title,tumor_content_estimate,
                     bam_tumor):
    '''
    '''
    args = [samtools_bin,"view",bam_tumor,"|", "head", "-n", "1", "|", "cut", "-f", "10"]                
    args=",".join(args).replace(',',' ')
    p1 =subprocess.Popen(args,stdout=subprocess.PIPE,shell=True)
    t = p1.communicate()[0]
    tumor_read_length=len(t.strip('\n'))
    print "tumor_read_length=",tumor_read_length
    if tumor_read_length == 0:
        logging.error("read length determination failed")
        return config.JOB_ERROR

    args = [Rscrip_dir,os.path.join(_pipeline_dir, "ExomeCNV.r"),
                    benign_cov,
                    tumor_cov,
                    loh_file,
                    cnv_file,
                    cnv_plot,
                    r_backup_file,
                    plot_title,
                    str(tumor_read_length),
                    tumor_content_estimate,
                    str(1)]#,"TRUE"]
    print args
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("Exome CNV calling failed")
        return config.JOB_ERROR
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("samtools_bin")
    parser.add_argument("Rscrip_dir")
    parser.add_argument("normal_cov")
    parser.add_argument("tumor_cov")
    parser.add_argument("loh_file")
    parser.add_argument("cnv_file")
    parser.add_argument("cnv_plot")
    parser.add_argument("r_backup_file")
    parser.add_argument("plot_title")
    parser.add_argument("tumor_content_estimate")
    parser.add_argument("bam_tumor")
    
    args = parser.parse_args()
    return run_exomeCNV(args.samtools_bin,args.Rscrip_dir, args.normal_cov, args.tumor_cov, 
                     args.loh_file, args.cnv_file,
                     args.cnv_plot,args.r_backup_file,
                     args.plot_title,args.tumor_content_estimate,
                     args.bam_tumor)


if __name__ == '__main__': 
    sys.exit(main())