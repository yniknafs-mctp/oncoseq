'''
Created on May 15, 2013

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import glob

import oncoseq.rnaseq.lib.config as config
from oncoseq.rnaseq.lib.libtable import Library

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="library_path_file", default=None)
    parser.add_argument("-d", "--dir", dest="input_dirs", action="append", default=None)
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
    if args.input_dirs is not None:
        for input_dir in args.input_dirs:
            if not os.path.exists(input_dir):
                parser.error("input dir '%s' not found" % (input_dir))
            if not os.path.isdir(input_dir):
                parser.error("input dir '%s' not a valid directory" % (input_dir))
            for path in glob.iglob(os.path.join(input_dir, "*")):
                if os.path.exists(path) and os.path.isdir(path):
                    library_paths.add(path)
    if len(library_paths) == 0:
        parser.error("No valid library paths found. Use '-i' to specify a file with "
                     "a list of result paths or '--dir' to search a directory for "
                     "results")     
    # read gene fusion information
    print '\t'.join(["sample_id", "library_id", "gtf_file", "bam_file"])
    for output_dir in sorted(library_paths):
        logging.debug("Processing %s" % (output_dir))
        library_xml_file = os.path.join(output_dir, config.LIBRARY_XML_FILE)
        if not os.path.exists(library_xml_file):
            logging.error("Library '%s' xml file not found" % (output_dir))
            continue
        library = list(Library.from_xml_file(library_xml_file))[0]
        config_xml_file = os.path.join(output_dir, config.CONFIG_XML_FILE)
        if not os.path.exists(config_xml_file):
            logging.error("Library '%s' config xml file not found" % (output_dir))
            continue
        # get results
        pipeline = config.PipelineConfig.from_xml(config_xml_file)
        results = config.RnaseqResults(library, pipeline, output_dir)
        if (not os.path.exists(results.tophat_bam_file) or 
            not os.path.exists(results.cufflinks_ab_initio_gtf_file)):
            logging.error("Library '%s' result file(s) not found" % (output_dir))
        else:
            fields = [library.sample_id, library.library_id, 
                      results.cufflinks_ab_initio_gtf_file,
                      results.tophat_bam_file]
            print '\t'.join(fields) 

if __name__ == '__main__':
    sys.exit(main())
