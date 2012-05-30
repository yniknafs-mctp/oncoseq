'''
Created on Jan 23, 2012

@author: mkiyer
'''
import argparse
import logging

# local imports
import pysam

# project imports
from oncoseq.lib.sam import parse_reads_by_qname

def filter_unmapped_pairs(bam_file,
                          r1_input_file, 
                          r2_input_file=None):
    input_sam_fhs = [pysam.Samfile(r1_input_file, "r")]
    bamfh = pysam.Samfile(bam_file, "wb", template=input_sam_fhs[0])
    if r2_input_file is not None:
        input_sam_fhs.append(pysam.Samfile(r2_input_file, "r"))
    read_iters = [parse_reads_by_qname(f) for f in input_sam_fhs]
    try:
        while True:            
            pe_reads = [it.next() for it in read_iters]
            # if either read maps to the abundant sequences, discard the pair
            pe_mapped_reads = []
            for readnum, reads in enumerate(pe_reads):
                mapped_reads = []
                for r in reads:
                    if not r.is_unmapped:
                        r.is_paired = True
                        r.mate_is_unmapped = True
                        if readnum == 1:
                            r.is_read2 = True
                        else:
                            r.is_read1 = True
                        mapped_reads.append(r)
                pe_mapped_reads.append(mapped_reads)
            if not all((len(mapped_reads) == 0) for mapped_reads in pe_mapped_reads):
                for readnum,mapped_reads in enumerate(pe_mapped_reads):
                    if len(mapped_reads) == 0:
                        bamfh.write(pe_reads[readnum][0])
                    else:
                        for r in mapped_reads:
                            bamfh.write(r)
    except StopIteration:
        pass

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file")
    parser.add_argument("read1_input_file")
    parser.add_argument("read2_input_file", nargs="?", default=None)
    args = parser.parse_args()
    filter_unmapped_pairs(args.bam_file,
                          args.read1_input_file,
                          args.read2_input_file)
    return 0
    
if __name__ == '__main__':
    main()