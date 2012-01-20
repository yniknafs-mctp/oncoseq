'''
Created on Aug 15, 2011

@author: mkiyer
'''
import argparse
import logging
import sys

from oncoseq.lib import config

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("output_file")
    args = parser.parse_args()
    open(args.output_file, "w").close()
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())