'''
Created on Dec 2, 2012

@author: mkiyer
'''
import logging
import argparse
import os
import sys
import subprocess

import pysam

from oncoseq.rnaseq.lib.base import parse_sam, remove_multihits, \
    to_fastq, check_executable
import oncoseq.rnaseq.pipeline
_pipeline_dir = oncoseq.rnaseq.pipeline.__path__[0]

def bam_to_fastq(bam_file, fastq_prefix, assume_sorted, 
                 tmp_dir): 
    # try to get the sort order from the header
    sort_order = None
    try:
        f = pysam.Samfile(bam_file, 'rb')
        if 'HD' in f.header:
            hd_dict = f.header['HD']
            if 'SO' in hd_dict:
                sort_order = hd_dict['SO']        
        f.close()
    except Exception as e:
        logging.warning("Error while trying to read BAM header: %s" % (str(e)))
        sort_order = None
    # decide whether to sort bam file
    if assume_sorted or (sort_order == "queryname"):
        sorted_bam_file = bam_file
        keep_sorted_bam = True
    else:
        # sort bam by queryname
        prefix = os.path.join(tmp_dir, os.path.splitext(os.path.basename(bam_file))[0])
        sorted_bam_prefix = prefix + '.qnamesorted'
        sorted_bam_file = prefix + '.qnamesorted.bam'
        keep_sorted_bam = False
        logging.info("sorting input file by queryname")
        args = ["samtools", "sort", "-m", "2000000000", 
                "-n", bam_file, sorted_bam_prefix]
        logging.debug("bam sort args: %s" % (map(str, args)))    
        retcode = subprocess.call(args)
        if retcode != 0:
            logging.error("bam sort failed")
            if os.path.exists(sorted_bam_file):
                os.remove(sorted_bam_file)
            return 1
        logging.debug("bam sort done")
    # extract paired and unpaired reads
    # setup output files
    paired_files = []
    paired_fhs = []
    paired_counts = 0
    unpaired_files = []
    unpaired_fhs = []
    unpaired_counts = [0,0]
    for readnum in (1,2):
        filename = "%s.paired.%d.fq" % (fastq_prefix, readnum)
        paired_files.append(filename)
        paired_fhs.append(open(filename, 'w'))
        filename = "%s.unpaired.%d.fq" % (fastq_prefix, readnum)
        unpaired_files.append(filename)
        unpaired_fhs.append(open(filename, 'w'))
    # parse bam
    logging.debug("writing fastq files")
    bamfh = pysam.Samfile(sorted_bam_file, 'rb')
    for pe_reads in parse_sam(bamfh, remove_suffix=True):
        r1, r2 = remove_multihits(pe_reads)
        if (r1 is None) and (r2 is None):
            raise Exception("Error parsing SAM file")
        elif (r1 is not None) and (r2 is not None):
            # pair
            print >>paired_fhs[0], to_fastq(r1, 0)
            print >>paired_fhs[1], to_fastq(r2, 1)
            paired_counts += 1
        else:
            if (r1 is not None):
                # unpaired read 1
                print >>unpaired_fhs[0], to_fastq(r1, 0)
                unpaired_counts[0] += 1
            elif (r2 is not None):
                # unpaired read 2
                print >>unpaired_fhs[1], to_fastq(r2, 1)
                unpaired_counts[1] += 1
    bamfh.close()
    logging.debug("found %d paired reads" % (paired_counts))
    logging.debug("found %d unpaired read1" % (unpaired_counts[0]))
    logging.debug("found %d unpaired read2" % (unpaired_counts[1]))
    # clean up
    for readnum in xrange(2):
        paired_fhs[readnum].close()
        unpaired_fhs[readnum].close()
    if not keep_sorted_bam:
        if os.path.exists(sorted_bam_file):
            os.remove(sorted_bam_file)
    return 0

def main():
    # setup logging
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # command line parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--assume-sorted', dest="assume_sorted",
                        action="store_true", default=False)
    parser.add_argument('--tmp-dir', dest="tmp_dir", default="/tmp")
    parser.add_argument('bam_file')
    parser.add_argument('fastq_prefix')
    args = parser.parse_args()
    if not os.path.exists(args.bam_file):
        parser.error("Input bam file not found")
    if not os.path.exists(args.tmp_dir):
        parser.error("tmp dir %s not found" % (args.tmp_dir))
    # check command line args
    if not check_executable("samtools"):
        parser.error("samtools binary not found")
    tmp_dir = os.path.abspath(args.tmp_dir)
    return bam_to_fastq(args.bam_file, 
                        args.fastq_prefix, 
                        args.assume_sorted,
                        tmp_dir)

if __name__ == '__main__':
    sys.exit(main())
