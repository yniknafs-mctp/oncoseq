'''
Created on Aug 9, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import sys

from oncoseq.lib import config
from oncoseq.lib import fastqc
from oncoseq.lib.base import file_exists_and_nz_size
from oncoseq.lib.seq import get_qual_conversion_func, FASTQ_QUAL_FORMATS, open_compressed

def parse_lines(line_iter, numlines=1):
    """
    generator that returns list of 'numlines' lines at a time
    """
    try:
        while True:
            yield [line_iter.next().rstrip() for x in xrange(numlines)]
    except StopIteration:
        pass

def copy_uncompress_sequence(src, dst, quals, fastqc_data_files):
    """
    uncompresses reads, renames reads, and converts quality scores 
    to 'sanger' format
    """
    # get fastq format from fastqc output
    encodings = set()
    for f in fastqc_data_files:
        if file_exists_and_nz_size(f):
            encoding = fastqc.get_fastq_encoding(f)
            if encoding not in fastqc.ENCODING_TO_QUAL_FORMAT:
                logging.error("Unrecognized FASTQ encoding %s" % (encoding))
                return config.JOB_ERROR
            encodings.add(fastqc.ENCODING_TO_QUAL_FORMAT[encoding])
    if len(encodings) > 1:
        logging.error("Detected different FASTQ encodings in paired-end files")
        return config.JOB_ERROR
    elif len(encodings) == 0:
        logging.warning("Could not locate FASTQC data to determine encoding, using value %s found in XML" % (quals))
        encodings.add(quals)
    # reconcile disagreement in quality scores
    auto_quals = encodings.pop()
    if quals != auto_quals:
        logging.warning("Auto-detected FASTQ format %s differs from XML file (%s)"
                        % (auto_quals, quals))
        quals = auto_quals
    logging.info("Copying and converting FASTQ files with quality scores "
                 "in %s format" % (quals))
    # setup file iterators
    fh = open_compressed(src)
    fqiter = parse_lines(fh, numlines=4)
    outfh = open(dst, "w")
    qual_func = get_qual_conversion_func(quals)
    linenum = 0    
    try:
        while True:
            lines = fqiter.next()
            # ignore extra tag name information
            lines[0] = lines[0].split()[0]
            # ignore redundant header
            lines[2] = "+"
            # convert quality score to sanger
            lines[3] = qual_func(lines[3])
            print >>outfh, '\n'.join(lines)
            linenum += 1
    except StopIteration:
        pass
    except:
        logging.error("Unexpected error during FASTQ file processing")
        if os.path.exists(dst):
            os.remove(dst)
        return config.JOB_ERROR
    fh.close()
    logging.debug("Processed %d fragments" % (linenum))
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--quals", choices=FASTQ_QUAL_FORMATS)
    parser.add_argument("--fastqc-data-files", dest="fastqc_data_files", default=None)
    parser.add_argument("src")
    parser.add_argument("dst")
    args = parser.parse_args()
    if args.fastqc_data_files is None:
        fastqc_data_files = []
    else:
        fastqc_data_files = args.fastqc_data_files.split(",")
    return copy_uncompress_sequence(args.src, args.dst, args.quals, fastqc_data_files)

if __name__ == '__main__': 
    sys.exit(main())