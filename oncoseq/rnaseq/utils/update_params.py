'''
Created on Apr 29, 2013

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import glob

import subprocess
import shutil

import oncoseq.rnaseq.lib.config as config
from oncoseq.rnaseq.lib.base import detect_format
from oncoseq.rnaseq.lib.libtable import Library, read_library_table_xls, FRAGMENT_LAYOUT_PAIRED
import oncoseq.rnaseq.pipeline
_pipeline_dir = oncoseq.rnaseq.pipeline.__path__[0]


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="library_path_file", default=None)
    parser.add_argument("-d", "--dir", dest="input_dirs", action="append", default=None)
    parser.add_argument("--backup", dest="backup", action="store_true", default=False)
    inp_group = parser.add_mutually_exclusive_group(required=True)
    inp_group.add_argument("--xls", dest="library_xls_file", 
                           default=None,
                           help="Excel (.xls/.xlsx) spreadsheet "
                           "containing RNA-Seq library information")
    inp_group.add_argument("--xml", dest="library_xml_file", 
                           default=None,
                           help="XML formatted file containing RNA-Seq "
                           "library information")
    args = parser.parse_args()
    # get library paths either from a file or by searching directories    
    library_paths = set()
    if args.library_path_file is not None:
        if not os.path.exists(args.library_path_file):
            parser.error("'%s' not found" % (args.library_path_file))
        for line in open(args.library_path_file):
            path = line.strip()
            if os.path.exists(path) and os.path.isdir(path):
                library_paths.add(path)
            else:
                logging.error("Directory not found at path %s" % (path))
    elif args.input_dirs is not None:
        for input_dir in args.input_dirs:
            if not os.path.exists(input_dir):
                parser.error("input dir '%s' not found" % (input_dir))
            if not os.path.isdir(input_dir):
                parser.error("input dir '%s' not a valid directory" % (input_dir))
            for path in glob.iglob(os.path.join(input_dir, "*")):
                if os.path.exists(path) and os.path.isdir(path):
                    library_paths.add(path)
    else:
        parser.error("Use '-i' to specify a file with a list of result "
                     "paths or '--dir' to search a directory for results")     
    # read libraries
    if args.library_xml_file is not None:
        if not os.path.exists(args.library_xml_file):
            parser.error("XML file %s not found" % (args.library_xml_file))
        # read library from XML
        new_libraries = {}
        for library in Library.from_xml_file(args.library_xml_file):
            new_libraries[library.library_id] = library
    elif args.library_xls_file is not None:                
        if not os.path.exists(args.library_xls_file):
            parser.error("Excel file %s not found" % (args.library_xls_file))
        # read libraries from XLS/XLSX
        logging.info("Reading library table '%s'" % (args.library_xls_file))
        new_libraries = read_library_table_xls(args.library_xls_file)
    # update library information
    for output_dir in library_paths:
        # read existing library information
        logging.debug("Processing %s" % (output_dir))
        library_xml_file = os.path.join(output_dir, config.LIBRARY_XML_FILE)
        library_id = os.path.basename((output_dir))
        # see if updated information exists
        if library_id in new_libraries:
            logging.info("Replacing library XML %s at %s" % (library_id, output_dir))
            new_libraries[library_id].to_xml_file(library_xml_file)
    
if __name__ == '__main__':
    sys.exit(main())
    
