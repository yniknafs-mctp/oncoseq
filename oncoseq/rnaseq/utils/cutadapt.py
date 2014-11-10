'''
Created on Nov 7, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import subprocess


CUTADAPT_BIN = '/mctp/projects/rnaseq/version_001_2013-01-14/sw/cutadapt/cutadapt-1.6/bin/cutadapt'


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    parser.add_argument("adapter1")
    parser.add_argument("adapter2")
    parser.add_argument("fastq1")
    parser.add_argument("fastq2")
    parser.add_argument("-m", dest='min', default = '50')
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.debug('Trimming Forward Reads')
    TMP1 = args.name + '.TMP1.fastq'
    TMP2 = args.name + '.TMP2.fastq'
    
    cut_args = [
                CUTADAPT_BIN, 
                '-m', args.min,
                '-a', args.adapter1,
                '-o', TMP1,
                '-p', TMP2, 
                args.fastq1,
                args.fastq2
                ]
    subprocess.call(cut_args)
    
    OUT1 = args.name + '_1.fastq'
    OUT2 = args.name + '_2.fastq'
    
    logging.debug('Trimming Reverse Reads')
    cut_args = [
                CUTADAPT_BIN, 
                '-m', args.min,
                '-a', args.adapter1,
                '-o', OUT2,
                '-p', OUT1, 
                TMP2,
                TMP1
                ]
    subprocess.call(cut_args)
    
    logging.debug('Removing temp files')
    os.remove(TMP1)
    os.remove(TMP2)
    
    logging.debug('Done')
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

