'''
Created on Oct 15, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import glob

import oncoseq.rnaseq.lib.config as config
import oncoseq.rnaseq.lib.picard as picard
import oncoseq.rnaseq.lib.fastqc as fastqc
from oncoseq.rnaseq.lib.libtable import Library

RESULT_HEADER_FIELDS = ["results_valid",
                        "total_fragments",
                        "total_aligned_reads", 
                        "mean_read_length", 
                        'aligned_bases',
                        'pct_mrna_bases',
                        'pct_intronic_bases',
                        'pct_intergenic_bases',
                        'median_5prime_to_3prime_bias',
                        'median_cv_coverage',
                        'bam_file', 
                        'gtf_file']

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="library_path_file", default=None)
    parser.add_argument("-d", "--dir", dest="input_dirs", action="append", default=None)
    args = parser.parse_args()
    # get library paths either from a file or by searching directories    
    library_paths = []
    if args.library_path_file is not None:
        if not os.path.exists(args.library_path_file):
            parser.error("'%s' not found" % (args.library_path_file))
        for line in open(args.library_path_file):
            path = line.strip()
            if os.path.exists(path) and os.path.isdir(path):
                library_paths.append(path)
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
                    library_paths.append(path)
    else:
        parser.error("Use '-i' to specify a file with a list of result "
                     "paths or '--dir' to search a directory for results")     
    # read libraries
    libraries = {}
    for output_dir in library_paths:
        library_xml_file = os.path.join(output_dir, config.LIBRARY_XML_FILE)
        if not os.path.exists(library_xml_file):
            logging.error("Library '%s' xml file not found" % (output_dir))
            continue
        library = list(Library.from_xml_file(library_xml_file))[0]
        #if not library.is_valid():
        #    logging.error("Library '%s' not valid" % (output_dir))
        #    continue
        # look for results
        config_xml_file = os.path.join(output_dir, config.CONFIG_XML_FILE)
        is_valid = False
        results = None
        if os.path.exists(config_xml_file):
            pipeline = config.PipelineConfig.from_xml(config_xml_file)
            results = config.RnaseqResults(library, pipeline, output_dir)
            is_valid, missing_files = results.validate()
        libraries[library.library_id] = (library, results, is_valid)
    # build list of all parameters
    params = set()
    for library in libraries.itervalues():
        params.update(library.params.keys())
    sorted_params = sorted(params)
    # output table
    header_fields = []
    header_fields.extend(Library.fields)
    header_fields.extend(sorted_params)
    header_fields.extend(RESULT_HEADER_FIELDS)
    print '\t'.join(header_fields)
    for library_id in sorted(libraries):
        library, results, results_valid = libraries[library_id]
        # convert lists to strings
        library.read1_files = ','.join(library.read1_files)
        library.read2_files = ','.join(library.read2_files)
        library.bam_files = ','.join(library.bam_files)
        # add library fields
        fields = []
        for field_name in Library.fields:
            fields.append(getattr(library, field_name))
        for param in sorted_params:
            fields.append(library.params.get(param, "na"))
        # add result fields
        fields.append('yes' if results_valid else 'no')
        if results is None:
            fields.extend(['na'] * len(RESULT_HEADER_FIELDS))
        else:
            if os.path.exists(results.fastqc_data_files[0]):
                total_frags = fastqc.get_total_sequences(results.fastqc_data_files[0])
                fields.append(total_frags)
            else:
                fields.append('na')
            if os.path.exists(results.alignment_summary_metrics):
                obj = picard.AlignmentSummaryMetrics(results.alignment_summary_metrics)
                fields.append(obj.get_total_reads())
                fields.append(obj.get_mean_read_length()) 
            else:
                fields.extend(['na', 'na'])
            if os.path.exists(results.rnaseq_metrics):
                metrics_dict = picard.get_rnaseq_metrics(results.rnaseq_metrics)
                fields.extend([metrics_dict['PF_ALIGNED_BASES'],
                               metrics_dict['PCT_MRNA_BASES'],
                               metrics_dict['PCT_INTRONIC_BASES'],
                               metrics_dict['PCT_INTERGENIC_BASES'],
                               metrics_dict['MEDIAN_5PRIME_TO_3PRIME_BIAS'],
                               metrics_dict['MEDIAN_CV_COVERAGE']])
            else:
                fields.extend(['na'] * 6)
            if os.path.exists(results.tophat_bam_file):
                fields.append(results.tophat_bam_file)
            else:
                fields.append('na')
            if os.path.exists(results.cufflinks_gtf_file):
                fields.append(results.cufflinks_gtf_file)
            else:
                fields.append('na')
        print '\t'.join(map(str, fields))
    return config.JOB_SUCCESS


    parser.add_argument("library_xls_file")
    parser.add_argument("root_dir")
    args = parser.parse_args()
    # read library file
    logging.info("Reading library table '%s'" % (args.library_xls_file))
    libraries = read_library_table_xls(args.library_xls_file)
    # build list of all parameters
    params = set()
    for library in libraries.itervalues():
        params.update(library.params.keys())
    sorted_params = sorted(params)
    # output table
    logging.info("Generating library table")
    header_fields = []
    header_fields.extend(Library.fields)
    header_fields.extend(sorted_params)
    header_fields.extend(RESULT_HEADER_FIELDS)
    print '\t'.join(header_fields)
    for library_id in sorted(libraries):
        library = libraries[library_id]
        library.read1_files = ','.join(library.read1_files)
        library.read2_files = ','.join(library.read2_files)
        library.bam_files = ','.join(library.bam_files)
        # look for results
        output_dir = os.path.join(args.root_dir, library.library_id)
        config_xml_file = os.path.join(output_dir, config.CONFIG_XML_FILE)
        is_valid = False
        results = None
        if not os.path.exists(config_xml_file):
            setattr(library, 'description', 'invalid')
        else:
            pipeline = config.PipelineConfig.from_xml(config_xml_file)
            results = config.RnaseqResults(library, pipeline, output_dir)
            is_valid, missing_files = results.validate()
            if not is_valid:
                setattr(library, 'description', 'invalid')
            else:
                setattr(library, 'description', 'valid')
        # add library fields
        fields = []
        for field_name in Library.fields:
            fields.append(getattr(library, field_name))
        for param in sorted_params:
            fields.append(library.params.get(param, "na"))
        # add result fields
        if results is None:
            fields.extend(['na'] * len(RESULT_HEADER_FIELDS))
        else:
            if os.path.exists(results.fastqc_data_files[0]):
                total_frags = fastqc.get_total_sequences(results.fastqc_data_files[0])
                fields.append(total_frags)
            else:
                fields.append('na')
            if os.path.exists(results.alignment_summary_metrics):
                obj = picard.AlignmentSummaryMetrics(results.alignment_summary_metrics)
                fields.append(obj.get_total_reads())
                fields.append(obj.get_mean_read_length()) 
            else:
                fields.extend(['na', 'na'])
            if os.path.exists(results.rnaseq_metrics):
                metrics_dict = picard.get_rnaseq_metrics(results.rnaseq_metrics)
                fields.extend([metrics_dict['PF_ALIGNED_BASES'],
                               metrics_dict['PCT_MRNA_BASES'],
                               metrics_dict['PCT_INTRONIC_BASES'],
                               metrics_dict['PCT_INTERGENIC_BASES'],
                               metrics_dict['MEDIAN_5PRIME_TO_3PRIME_BIAS'],
                               metrics_dict['MEDIAN_CV_COVERAGE']])
            else:
                fields.extend(['na'] * 6)
            if os.path.exists(results.tophat_bam_file):
                fields.append(results.tophat_bam_file)
            else:
                fields.append('na')
            if os.path.exists(results.cufflinks_gtf_file):
                fields.append(results.cufflinks_gtf_file)
            else:
                fields.append('na')
        print '\t'.join(map(str, fields))
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())