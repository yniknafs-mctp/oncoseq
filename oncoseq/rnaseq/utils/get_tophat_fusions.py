'''
Created on Apr 24, 2013

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
        #param_field = ','.join(['%s=%s' % (k,library.params[k]) 
        #                        for k in sorted(library.params.keys())])
        if 'cancer_progression' in library.params:
            progression_field = library.params['cancer_progression']
        else:
            progression_field = 'na'
        fields = [library.study_id, library.cohort_id, library.patient_id,
                  library.sample_id, progression_field]
        # get fusions
        if not os.path.exists(results.tophat_fusion_post_result_file):
            logging.error("Library '%s' fusion result file not found" % (output_dir))
        else:
            fusion_lines = [x.strip().split('\t') for x in open(results.tophat_fusion_post_result_file)]
            num_fusions = len(fusion_lines)
        if num_fusions == 0:
            fusion_lines = [[library.library_id] + ['na'] * 10]
        # make output fields
        fields.append(str(num_fusions))
        for fusion_fields in fusion_lines:
            print '\t'.join(fields + fusion_fields) 
    
if __name__ == '__main__':
    sys.exit(main())
