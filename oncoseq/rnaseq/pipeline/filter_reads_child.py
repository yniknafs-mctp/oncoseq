'''
Created on Jan 10, 2013

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import array

# local imports
import pysam

# project imports
from oncoseq.rnaseq.lib.base import parse_sr_reads, parse_sam, \
    remove_multihits, to_fastq, to_fastq_sr

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file")
    parser.add_argument("bam_file")
    parser.add_argument("counts_file")
    parser.add_argument("fastq_files")
    args = parser.parse_args()
    # open input files
    suffix = os.path.splitext(args.sam_file)[-1]
    mode = 'rb' if suffix == '.bam' else 'r'
    insamfh = pysam.Samfile(args.sam_file, mode)
    nrefs = len(insamfh.references)
    # open output files
    outbamfh = pysam.Samfile(args.bam_file, "wb", template=insamfh)
    fastq_files = args.fastq_files.split(',')
    fastq_fhs = [open(f,'w') for f in fastq_files]
    # init counts
    def make_array(n):
        return array.array('L', (0 for x in xrange(n)))
    counts = {'UU': make_array(nrefs),
              'CP': make_array(nrefs),
              'DP': make_array(nrefs),
              'UP': make_array(nrefs)}
    unmapped = 0
    error = False
    # process single and paired-end differently    
    if len(fastq_files) == 1:
        for reads in parse_sr_reads(insamfh):
            if len(reads) == 0:
                error = True
                break
            r = reads[0]
            if r.is_unmapped:
                print >>fastq_fhs[0], to_fastq_sr(r)
                unmapped += 1
            else:
                altype = r.opt('YT')
                counts[altype][r.tid] += 1            
                outbamfh.write(r)
    else:
        for pe_reads in parse_sam(insamfh, 
                                  readnum_in_qname=False, 
                                  remove_suffix=True):
            r1, r2 = remove_multihits(pe_reads)
            # check for unmatched reads or parsing errors
            if (r1 is None) or (r2 is None):
                error = True
                break
            # if both reads are unmapped then write to fastq
            if (r1.is_unmapped and r2.is_unmapped):
                print >>fastq_fhs[0], to_fastq(r1, 0)
                print >>fastq_fhs[1], to_fastq(r2, 1)
                unmapped += 2
            else:
                for r in (r1, r2):
                    if not r.is_unmapped:
                        # 'YT' is bowtie2 SAM tag that gives PE alignment 
                        # information
                        altype = r.opt('YT')
                        counts[altype][r.tid] += 1            
                    outbamfh.write(r)
    # cleanup
    for f in fastq_fhs:
        f.close()
    outbamfh.close()
    if error:
        for i in xrange(len(fastq_files)):
            os.remove(fastq_files[i])
        os.remove(args.bam_file)
        logging.error("filter_reads_child quit unexpectedly")
        return 1
    # print count file
    countfileh = open(args.counts_file, "w")
    print >>countfileh, '\t'.join(["#rname", "length", "paired", 
                                   "paired_discordant", "paired_orphan", 
                                   "unpaired", "total"])
    for i in xrange(len(insamfh.references)):
        paired = counts['CP'][i]
        paired_discordant = counts['DP'][i]
        paired_orphan = counts['UP'][i]
        unpaired = counts['UU'][i]
        total = paired + paired_discordant + paired_orphan + unpaired
        fields = [insamfh.references[i],
                  insamfh.lengths[i],
                  paired,
                  paired_discordant,
                  paired_orphan,
                  unpaired,
                  total]
        print >>countfileh, '\t'.join(map(str, fields))
    print >>countfileh, '\t'.join(['*', '0', '0', '0', '0', '0', 
                                   str(unmapped)])
    countfileh.close()
    insamfh.close()
    return 0

if __name__ == '__main__':
    sys.exit(main())