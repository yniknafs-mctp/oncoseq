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
from oncoseq.rnaseq.lib.inspect import RnaseqLibraryMetrics

RESULT_HEADER_FIELDS = ["results_valid",
                        "total_fragments",
                        "strand_fraction",
                        "predicted_library_type",
                        "frag_size_mean",
                        "frag_size_stdev",
                        "total_aligned_reads", 
                        "mean_read_length", 
                        'aligned_bases',
                        'pct_mrna_bases',
                        'pct_intronic_bases',
                        'pct_intergenic_bases',
                        'median_5prime_to_3prime_bias',
                        'median_cv_coverage',
                        'bam_file',
                        'ab_initio_gtf_file']

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
        logging.debug("Processing %s" % (output_dir))
        library_xml_file = os.path.join(output_dir, config.LIBRARY_XML_FILE)
        if not os.path.exists(library_xml_file):
            logging.error("Library '%s' xml file not found" % (output_dir))
            continue
        library = list(Library.from_xml_file(library_xml_file))[0]
        config_xml_file = os.path.join(output_dir, config.CONFIG_XML_FILE)
        fields = []
        if not os.path.exists(config_xml_file):
            fields.extend(['na'] * len(RESULT_HEADER_FIELDS))        
        else:
            # get results
            pipeline = config.PipelineConfig.from_xml(config_xml_file)
            results = config.RnaseqResults(library, pipeline, output_dir)
            # validate results
            results_valid, missing_files = results.validate()
            fields.append('yes' if results_valid else 'no')
            # add results
            if os.path.exists(results.fastqc_data_files[0]):
                total_frags = fastqc.get_total_sequences(results.fastqc_data_files[0])
                fields.append(total_frags)
            else:
                fields.append('na')
            if os.path.exists(results.library_metrics_file):
                obj = RnaseqLibraryMetrics.from_file(results.library_metrics_file)
                strand_frac = obj.read1_strand_fraction()
                predicted_library_type = obj.predict_library_type()
                frag_size_mean = int(round(obj.tlen_at_percentile(50.0),0))
                frag_size_stdev = int(round(obj.std(),0))
                fields.extend([strand_frac,
                               predicted_library_type,
                               frag_size_mean,
                               frag_size_stdev])
            else:
                fields.extend(['na'] * 4)
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
            if os.path.exists(results.cufflinks_ab_initio_gtf_file):
                fields.append(results.cufflinks_ab_initio_gtf_file)
            else:
                fields.append('na')
        libraries[library.library_id] = (library, fields)
    # build list of all parameters
    params = set()
    for library, result_fields in libraries.itervalues():
        params.update(library.params.keys())
    sorted_params = sorted(params)
    # output table
    header_fields = []
    header_fields.extend(Library.fields)
    header_fields.extend(sorted_params)
    header_fields.extend(RESULT_HEADER_FIELDS)
    print '\t'.join(header_fields)
    for library_id in sorted(libraries):
        library, result_fields  = libraries[library_id]
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
        fields.extend(result_fields)
        print '\t'.join(map(str, fields))
    return config.JOB_SUCCESS
    
if __name__ == '__main__':
    sys.exit(main())