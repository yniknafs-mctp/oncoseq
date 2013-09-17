'''
Created on Aug 21, 2013

@author: mkiyer
'''
import argparse
import logging
import sys
import os

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", dest="species", default="human")
    parser.add_argument("--library-type", dest="library_type", default=FR_UNSTRANDED)
    parser.add_argument("--library-protocol", dest="library_protocol", default=POLYA_TRANSCRIPTOME)
    parser.add_argument("--param", dest="param_list", action="append", default=None)
    parser.add_argument("--xml", dest="write_xml", action="store_true", default=False)
    parser.add_argument("--ignore", dest="ignore_file", default=None)
    parser.add_argument("cghub_xml_file")
    parser.add_argument("seq_repo")    
    args = parser.parse_args()
    
if __name__ == '__main__':
    sys.exit(main())   
