'''
Created on Aug 9, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import sys

from pantyhose.lib import config
from pantyhose.lib.seq import get_qual_conversion_func, FASTQ_QUAL_FORMATS, open_compressed

def parse_lines(line_iter, numlines=1):
    """
    generator that returns list of 'numlines' lines at a time
    """
    try:
        while True:
            yield [line_iter.next().rstrip() for x in xrange(numlines)]
    except StopIteration:
        pass

def copy_uncompress_sequence(src, dst, quals):
    """
    uncompresses reads, renames reads, and converts quality scores 
    to 'sanger' format
    """
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
    parser.add_argument("src")
    parser.add_argument("dst")
    args = parser.parse_args()
    return copy_uncompress_sequence(args.src, args.dst, args.quals)
    
if __name__ == '__main__': 
    sys.exit(main())