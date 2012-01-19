'''
Created on Aug 14, 2011

@author: mkiyer
'''
import argparse
import logging
import sys

from oncoseq.lib import config
from oncoseq.lib.fragment_size_distribution import FragmentSizeDistribution

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("output_file")
    parser.add_argument("input_files", nargs="+")
    args = parser.parse_args()
    frag_size_dist = None
    for input_file in args.input_files:
        lane_frag_size_dist = FragmentSizeDistribution.from_file(open(input_file))
        if frag_size_dist is None:
            frag_size_dist = lane_frag_size_dist
        else:
            frag_size_dist = FragmentSizeDistribution.merge(frag_size_dist, lane_frag_size_dist)
    if frag_size_dist is None:
        logging.error("Library %d: no lanes available to provide fragment size distribution")
        return config.JOB_ERROR
    frag_size_dist.to_file(open(args.output_file, "w"))
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())