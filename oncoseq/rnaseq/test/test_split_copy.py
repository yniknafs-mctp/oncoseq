'''
Created on Jan 16, 2012

@author: mkiyer
'''
import logging
import sys
import argparse

from oncoseq.lib.cluster import remote_copy_file

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("src")
    parser.add_argument("dst")
    parser.add_argument("remote_address")
    parser.add_argument("username")
    parser.add_argument("port")
    parser.add_argument("tmp_dir")
    args = parser.parse_args()
    maxsize = (256<<20)
    print maxsize
    remote_copy_file(args.src, args.dst, args.remote_address, args.username, args.port, 
                     maxsize=maxsize, tmp_dir=args.tmp_dir)
    return 0

if __name__ == '__main__': 
    sys.exit(main())