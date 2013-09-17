'''
Created on Apr 30, 2013

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import collections
import shutil
import numpy as np

def parse_gtf_metadata_table(filename):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    tracking_id_ind = header_fields.index('tracking_id')
    tracking_ids = []
    for line in fileh:
        fields = line.strip().split('\t')
        tracking_ids.append(fields[tracking_id_ind])
    fileh.close()
    return tracking_ids

def parse_cufflinks_file(filename):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    tracking_id_ind = header_fields.index('tracking_id')
    fpkm_ind = header_fields.index('FPKM')
    data_dict = {}
    for line in fileh:
        fields = line.strip().split('\t')
        tracking_id = fields[tracking_id_ind]
        fpkm = float(fields[fpkm_ind])
        data_dict[tracking_id] = fpkm
    fileh.close()
    return data_dict

def parse_library_table(library_table_file, output_dir):
    outfile = os.path.join(output_dir, "library_table.txt")
    outfileh = open(outfile, 'w')
    fileh = open(library_table_file)
    header_fields = fileh.next().strip().split('\t')
    print >>outfileh, '\t'.join(header_fields)
    library_id_ind = header_fields.index('library_id')
    result_dir_ind = header_fields.index('results_dir')
    library_results = collections.OrderedDict()
    for line in fileh:
        if line.startswith("#"):
            continue
        if not line:
            continue
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        library_id = fields[library_id_ind]
        result_dir = fields[result_dir_ind]
        genes_file = os.path.join(result_dir, 'cufflinks_known', 'genes.fpkm_tracking')
        if not os.path.exists(genes_file):
            continue
        print >>outfileh, '\t'.join(fields)        
        library_results[library_id] = genes_file
    fileh.close()
    outfileh.close()
    return library_results

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_metadata_table")
    parser.add_argument("library_table")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # create output dir
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    # read gtf_metadata
    print 'reading gene metadata'
    feature_ids = parse_gtf_metadata_table(args.gtf_metadata_table)
    # copy to output directory
    shutil.copy(args.gtf_metadata_table, 
                os.path.join(args.output_dir, "gene_table.txt"))
    # read library table
    print 'reading library metadata'
    library_results = parse_library_table(args.library_table, args.output_dir)
    # create matrix
    print 'found %d genes and %d libraries' % (len(feature_ids), len(library_results))
    matrix_file = os.path.join(args.output_dir, 'gene_expression.mmap')
    mat = np.memmap(matrix_file, dtype='float32', mode='w+', 
                    shape=(len(feature_ids),len(library_results)))
    # fill matrix with data
    for j,library_id in enumerate(library_results.iterkeys()):
        print j, library_id
        genes_file = library_results[library_id]
        gene_exp_dict = parse_cufflinks_file(genes_file)
        for feature_id in feature_ids:
            if feature_id not in gene_exp_dict:
                logging.warning("Cannot find %s in library %s" % (feature_id, library_id))
        gene_exp_list = [gene_exp_dict.get(x, 0.0) for x in feature_ids]
        mat[:,j] = np.array(gene_exp_list, dtype='float32')
    # write to file
    fileh = open(os.path.join(args.output_dir, 'gene_expression.txt'), 'w')
    header_fields = ['feature_id']
    header_fields.extend(library_results.keys())
    print >>fileh, '\t'.join(header_fields)
    for i,feature_id in enumerate(feature_ids):
        fields = [feature_id]
        fields.extend(map(str, mat[i,:]))
        print >>fileh, '\t'.join(fields)
    fileh.close()    
    return 0

if __name__ == '__main__': 
    sys.exit(main())