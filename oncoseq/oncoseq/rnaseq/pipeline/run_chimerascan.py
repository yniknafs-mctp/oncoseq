'''
Created on Feb 25, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import subprocess

from oncoseq.lib import config
from oncoseq.lib.fragment_size_distribution import FragmentSizeDistribution
from oncoseq.lib.seq import detect_read_length

def run_chimerascan(output_dir, index, fastq_files, 
                    frag_size_dist_file, frag_size_percentile, 
                    trim5, trim3, library_type, num_processors,
                    extra_args, chimerascan_bin):
    # get read length and fragment size parameters
    orig_read_length = detect_read_length(fastq_files[0])
    read_length = orig_read_length - trim5 - trim3
    frag_size_dist = FragmentSizeDistribution.from_file(open(frag_size_dist_file))
    # choose a segment length to optimize mapping
    optimal_frag_size = frag_size_dist.isize_at_percentile(frag_size_percentile)
    logging.info("Fragment size at %f percent of distribution is %d" % 
                 (frag_size_percentile, optimal_frag_size))
    segment_length = int(round(optimal_frag_size / 2.0))
    logging.info("Optimal segment length is %d" % (segment_length))
    segment_length = min(segment_length, read_length)
    segment_length = max(config.CHIMERASCAN_MIN_SEGMENT_LENGTH, segment_length)
    logging.info("After adjusting for min %d and read length %d, final segment length is %d" % 
                 (config.CHIMERASCAN_MIN_SEGMENT_LENGTH, read_length, segment_length))
    # setup run
    args = [chimerascan_bin,
            "-p", num_processors,
            "--quals", "sanger",
            "--library-type", library_type,
            "--trim5", trim5,
            "--trim3", trim3,
            "--segment-length", segment_length]
    for arg in extra_args:
        args.extend(arg.split())
    args.append(index)
    args.extend(fastq_files)
    args.append(output_dir)
    logging.debug("\targs: %s" % (map(str, args)))
    retcode = subprocess.call(map(str, args))
    return retcode

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", dest="bin", default="chimerascan_run.py")
    parser.add_argument("--arg", dest="extra_args", action="append", default=[])
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--library-type", dest="library_type", default="fr-unstranded")
    parser.add_argument("--trim5", type=int, dest="trim5")
    parser.add_argument("--trim3", type=int, dest="trim3")
    parser.add_argument("--frag-size-percentile", type=float, dest="frag_size_percentile")
    parser.add_argument("output_dir")
    parser.add_argument("index")
    parser.add_argument("frag_size_dist_file")
    parser.add_argument("fastq_files", nargs="+")
    args = parser.parse_args()
    run_chimerascan(output_dir=args.output_dir,
                    index=args.index,
                    fastq_files=args.fastq_files, 
                    frag_size_dist_file=args.frag_size_dist_file,
                    frag_size_percentile=args.frag_size_percentile, 
                    trim5=args.trim5, 
                    trim3=args.trim3,
                    library_type=args.library_type, 
                    num_processors=args.num_processors,
                    extra_args=args.extra_args,
                    chimerascan_bin=args.bin)
    
if __name__ == '__main__':
    sys.exit(main())