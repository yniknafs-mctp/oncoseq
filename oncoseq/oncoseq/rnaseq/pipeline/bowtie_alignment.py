'''
Created on Aug 9, 2011

@author: mkiyer
'''
import sys
import os
import logging
import subprocess
import argparse

from pantyhose.lib.seq import SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT
from pantyhose.lib import config

import pantyhose.pipeline
_pipeline_dir = pantyhose.pipeline.__path__[0]
_fix_order_script = os.path.join(_pipeline_dir, "fix_alignment_ordering.py")

translate_quals = {SOLEXA_FORMAT: 'solexa-quals',
                   ILLUMINA_FORMAT: 'solexa1.3-quals',
                   SANGER_FORMAT: 'phred33-quals'}

def align_sr(fastq_file, 
             bowtie_index,
             output_file, 
             unaligned_fastq_file=None,
             maxmultimap_fastq_file=None,
             trim5=0,
             trim3=0,
             num_processors=2, 
             quals=SANGER_FORMAT,
             maxhits=1,
             maxmulti=None,
             mismatches=2, 
             bowtie_bin="bowtie", 
             log_file=None,
             keep_unmapped=False,
             keep_maxmulti=False):
    if mismatches is None:
        mismatches = 2
    bowtie_threads = max(1, num_processors-1)
    args = [bowtie_bin, "-q", "-S", 
            "--%s" % translate_quals[quals],
            "-p", str(bowtie_threads),
            "-n", str(mismatches),
            "--trim5", trim5,
            "--trim3", trim3]    
    if (maxhits is None) and (maxmulti is None):
        args.append("-a")
    else:
        if maxhits is not None:
            args.extend(["-k", maxhits])
        if maxmulti is not None:
            args.extend(["-m", str(maxmulti)])
    if unaligned_fastq_file is not None:
        args.extend(["--un", unaligned_fastq_file])
    if maxmultimap_fastq_file is not None:
        args.extend(["--max", maxmultimap_fastq_file])    
    args.extend([bowtie_index, fastq_file])
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    # setup logging
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None
    # start bowtie alignment process
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # fix alignment ordering and convert to BAM, also extend sequences
    # back to full length by adding padding to CIGAR string
    args = [sys.executable, _fix_order_script]
    if maxmulti is not None:
        args.extend(["--maxmulti", str(maxmulti)])
    if keep_unmapped:
        args.append("--keep-unmapped")
    if keep_maxmulti:
        args.append("--keep-maxmulti")
    args.extend([fastq_file, "-", output_file])
    logging.debug("SAM alignment reordering args: %s" % (' '.join(args)))
    fix_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    # wait for processes to complete
    retcode1 = fix_p.wait()
    if retcode1 != 0:
        logging.error("SAM to BAM conversion script failed")
        # kill alignment process
        aln_p.kill()
        # cleanup output file
        if os.path.exists(output_file):
            os.remove(output_file)
        # end logging
        if logfh is not None:
            logfh.close()
        return config.JOB_ERROR
    retcode2 = aln_p.wait()
    # end logging
    if logfh is not None:
        logfh.close()
    if retcode2 != 0:
        logging.error("Alignment process failed")
        # cleanup output file
        if os.path.exists(output_file):
            os.remove(output_file)
        return config.JOB_ERROR
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=2)
    parser.add_argument("--un", dest="unaligned_fastq_file", default=None)
    parser.add_argument("--mmax", dest="maxmultimap_fastq_file", default=None)
    parser.add_argument("--trim5", type=int, dest="trim5", default=0)
    parser.add_argument("--trim3", type=int, dest="trim3", default=0)
    parser.add_argument("--quals", dest="quals")
    parser.add_argument("--maxhits", type=int, dest="maxhits", default=1)
    parser.add_argument("--maxmulti", type=int, dest="maxmulti", default=None)
    parser.add_argument("--mismatches", type=int, dest="mismatches", default=2)
    parser.add_argument("--bowtie-bin", dest="bowtie_bin", default="bowtie")
    parser.add_argument("--log", dest="log_file")
    parser.add_argument("--keep-unmapped", action="store_true", dest="keep_unmapped", default=False)
    parser.add_argument("--keep-maxmulti", action="store_true", dest="keep_maxmulti", default=False)
    parser.add_argument("fastq_file")
    parser.add_argument("bowtie_index")
    parser.add_argument("output_file")
    args = parser.parse_args()
    return align_sr(args.fastq_file, 
                    args.bowtie_index,
                    args.output_file, 
                    unaligned_fastq_file=args.unaligned_fastq_file,
                    maxmultimap_fastq_file=args.maxmultimap_fastq_file,
                    trim5=args.trim5,
                    trim3=args.trim3,
                    num_processors=args.num_processors, 
                    quals=args.quals,
                    maxhits=args.maxhits,
                    maxmulti=args.maxmulti,
                    mismatches=args.mismatches,
                    bowtie_bin=args.bowtie_bin,
                    log_file=args.log_file,
                    keep_unmapped=args.keep_unmapped,
                    keep_maxmulti=args.keep_maxmulti)
    
if __name__ == '__main__':
    main()
