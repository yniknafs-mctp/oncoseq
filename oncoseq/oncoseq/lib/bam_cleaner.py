'''
Created on Mar 2, 2012

@author: oabalbin
'''
import os
import sys
import logging
import tempfile
import subprocess
import argparse

from oncoseq.lib.base import up_to_date
from oncoseq.lib import config


def markDuplicates(indexed_bam_file, outfile, picard_dir, use_mem=2048):
    '''
    '''
    picard_command=os.path.join(picard_dir,'MarkDuplicates.jar')
    metrics_file=indexed_bam_file.replace('.bam','.duplicates.stats')
    
    args=['java','-Xmx'+str(use_mem)+'m','-jar', picard_command, 
          'I='+indexed_bam_file, 'O='+outfile, 'METRICS_FILE='+metrics_file,
          'ASSUME_SORTED=true','VALIDATION_STRINGENCY=SILENT',
          'REMOVE_DUPLICATES=true'
          ]
    
    cmd = ",".join(args).replace(',',' ').replace(';',',')
    return cmd


def sortIndexSam(input_sam_file, output_bam_file, 
                picard_dir, sort_order='coordinate',
                use_mem=2048):
    '''
    sort sam or bam
    SO=SortOrder Sort order of output file  Required. Possible values: {unsorted, queryname, coordinate} 
    '''
    picard_command = os.path.join(picard_dir,"SortSam.jar")

    args = ['java', '-Xmx'+str(use_mem)+'m', '-jar', picard_command, 
            'I=%s' % (input_sam_file), 
            'O=%s' % (output_bam_file),
            'SORT_ORDER=%s' % (sort_order),
            'CREATE_INDEX=true',
            'VALIDATION_STRINGENCY=SILENT']
    
    cmd = ",".join(args).replace(',',' ').replace(';',',')
    return cmd


def bam_cleaner(bam_file, bam_smdup_file, picard_dir,tmp_dir):
    '''
    It cleans and merged all the lane files belonging to the 
    same sample
    submit_job_func
    '''
    # Tmp Files
    odir = os.path.dirname(tmp_dir)
    fh, btmp = tempfile.mkstemp(suffix='bam', prefix='tmp', dir=odir)
    os.close(fh)
    # cluster parameters        
    # Lists
                        
    if up_to_date(bam_smdup_file,bam_smdup_file):
        logging.info("[SKIPPED] SAM BAM CLEANER. Mark duplicates step. File %s file is up to date" % (bam_smdup_file))
    else:       
        # Mark duplicates in the merged Bam file\n
        logging.info("SAM BAM CLEANER. Remove duplicates step for file %s" % (bam_smdup_file))
        args = markDuplicates(bam_file, btmp,picard_dir) #,MEM_PER_CORE,
        print args
        retcode = subprocess.call(args,shell=True)
        
        if retcode != 0:
            logging.error("SAM BAM CLEANER: Removing duplicates failed")
            if os.path.exists(btmp):
                os.remove(btmp)
            return config.JOB_ERROR
        
        args  = sortIndexSam(btmp, bam_smdup_file, picard_dir) #,MEM_PER_CORE)
        print args
        retcode = subprocess.call(args,shell=True)

        if retcode != 0:
            logging.error("SAM BAM CLEANER: Sorting BAM without duplicates file failed")
            if os.path.exists(btmp):
                os.remove(btmp)
            return config.JOB_ERROR
        
        args = ['rm',btmp]
        retcode = subprocess.call(args)

        if retcode != 0:
            logging.error("SAM BAM CLEANER: Removing temporary files failed")
            if os.path.exists(btmp):
                os.remove(btmp)
            return config.JOB_ERROR
            
    return  config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--picard-dir", dest="picard_dir", default="")
    parser.add_argument("--tmp-dir", dest="tmp_dir", default="")    
    parser.add_argument("bam_file")
    parser.add_argument("bam_smdup")
    args = parser.parse_args()
    
    return bam_cleaner(args.bam_file,
                       args.bam_smdup,
                       args.picard_dir, 
                       args.tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())
