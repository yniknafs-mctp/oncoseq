'''
Created on Jun 2, 2011

@author: mkiyer
'''
import logging
import collections
import os
import argparse

# local imports
import pysam

ReorderBufferItem = collections.namedtuple('ReorderBufferItem', ("fqrec", "reads"))
FASTQRecord = collections.namedtuple('FASTQRecord', ("qname", "seq", "qual"))

def parse_fastq_record(line_iter):
    try:
        qname = line_iter.next().rstrip()[1:].split()[0]
        seq = line_iter.next().rstrip()
        line_iter.next()
        qual = line_iter.next().rstrip()
        yield FASTQRecord(qname, seq, qual)
        while True:
            qname = line_iter.next().rstrip()[1:].split()[0]
            seq = line_iter.next().rstrip()
            line_iter.next()
            qual = line_iter.next().rstrip()
            yield FASTQRecord(qname, seq, qual)
    except StopIteration:
        pass

def fix_sr_alignment_ordering(samfh, fastqfh, maxlen=100000):
    # setup fastq file iterator 
    fqiter = parse_fastq_record(fastqfh)
    # function for initializing new buffer entry
    buf_init_func = lambda fqrec: ReorderBufferItem(fqrec, [])
    # initialize the qname dictionary to match the fastq file    
    buf = collections.deque()
    qname_read_dict = {}
    for read in samfh:
        # set key for indexing reads
        key = read.qname
        # check if this read is already in the buffer
        if key not in qname_read_dict:
            # if buffer full empty the first entries
            while len(buf) >= maxlen:
                # get first key in buf
                first_key = buf.popleft()
                # return reads at this qname, then delete them
                yield qname_read_dict[first_key]
                del qname_read_dict[first_key]
            # add new qnames to buffer
            while True:                
                # get next qname from fastq file and add it to the queue
                fqrec = fqiter.next()
                next_key = fqrec.qname
                buf.append(next_key)
                qname_read_dict[next_key] = buf_init_func(fqrec)
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_key == key:
                    break
        # add read to buffer
        qname_read_dict[key].reads.append(read)
    # empty remaining entries in buffer
    while len(buf) > 0:
        yield qname_read_dict[buf.popleft()]

def reorder_sam(input_fastq_file, input_file, output_file,
                maxmulti=None, 
                keep_unmapped=True,
                keep_maxmulti=True,
                input_bam=False,
                output_bam=False):
    # detect file formats and open input and output files
    input_ext = os.path.splitext(input_file)[1]
    output_ext = os.path.splitext(output_file)[1]
    if input_ext == ".sam":
        input_mode = "r"
    elif input_ext == ".bam":
        input_mode = "rb"
    else:
        input_mode = "rb" if input_bam else "r"
    if output_ext == ".sam":
        output_mode = "wh"
    elif output_ext == ".bam":
        output_mode = "wb"
    else:
        output_mode = "wb" if output_bam else "w"
    infh = pysam.Samfile(input_file, input_mode)
    fastqfh = open(input_fastq_file, "r")
    outfh = pysam.Samfile(output_file, output_mode, template=infh)
    reorder_func = fix_sr_alignment_ordering(infh, fastqfh)
    # iterate through buffer
    num_frags = 0
    num_unmapped = 0
    num_maxmulti = 0
    for bufitem in reorder_func:
        num_frags += 1
        for r in bufitem.reads:
            # keep statistics of unmapped/multimapped reads and
            # suppress output if 'keep_unmapped' is False
            if r.is_unmapped:
                xm_tag = r.opt('XM')
                if (maxmulti is None) or (xm_tag < maxmulti):
                    num_unmapped += 1
                    if not keep_unmapped:
                        continue
                else:
                    num_maxmulti += 1
                    if not keep_maxmulti:
                        continue
            outfh.write(r)
    outfh.close()
    infh.close()
    fastqfh.close()
    logging.debug("Found %d fragments" % (num_frags))
    logging.debug("\t%d unmapped reads" % (num_unmapped))
    logging.debug("\t%d multimapping (>=%s) reads" % 
                  (num_maxmulti, 
                   "inf" if maxmulti is None else str(maxmulti)))


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxmulti", type=int, dest="maxmulti", default=None)
    parser.add_argument("--keep-unmapped", action="store_true", dest="keep_unmapped", default=False)
    parser.add_argument("--keep-maxmulti", action="store_true", dest="keep_maxmulti", default=False)
    parser.add_argument("--input-bam", action="store_true", dest="input_bam", default=False)
    parser.add_argument("--output-bam", action="store_true", dest="output_bam", default=False)
    parser.add_argument("fastq_file")
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    reorder_sam(args.fastq_file, args.input_file, args.output_file,
                maxmulti=args.maxmulti, 
                keep_unmapped=args.keep_unmapped,
                keep_maxmulti=args.keep_maxmulti,
                input_bam=False,
                output_bam=False)
    
if __name__ == '__main__':
    main()