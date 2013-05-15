'''
Created on May 1, 2013

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import numpy as np

GENE_EXPRESSION_DIR = "/mctp/users/mkiyer/rnaseq/version_001_2013-01-14/data_freeze_2013-04-24/gene_expression"

from get_gene_expression_table import read_lines, read_table

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("library_ids")
    parser.add_argument("gene_ids")
    parser.add_argument("pheno_file")
    parser.add_argument("expr_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # read input library ids
    library_ids = read_lines(args.library_ids)
    library_id_set = set(library_ids)
    # read input gene ids
    gene_ids = read_lines(args.gene_ids)
    gene_id_set = set(gene_ids)
    # find libraries
    library_table_file = os.path.join(GENE_EXPRESSION_DIR, "library_table.txt")
    ncols, lib_header_fields, lib_metadata, lib_inds, lib_ids = \
        read_table(library_table_file, primary_key='library_id', 
                   subset=library_id_set)
    if len(lib_inds) == 0:
        print "No libraries found"
        return 1
    # find genes
    gene_table_file = os.path.join(GENE_EXPRESSION_DIR, "gene_table.txt")
    nrows, g_header_fields, g_metadata, g_inds, g_ids = \
        read_table(gene_table_file, primary_key='tracking_id', 
                   subset=gene_id_set)
    if len(g_inds) == 0:
        print "No genes found"
        return 1
    # read gene expression
    matrix_file = os.path.join(GENE_EXPRESSION_DIR, "gene_expression.mmap")
    mat = np.memmap(matrix_file, dtype='float32', mode='r', 
                    shape=(nrows,ncols))
    # get subset of matrix
    submat = mat[g_inds,:]
    submat = submat[:,lib_inds]
    # write pheno file
    fileh = open(args.pheno_file, 'w')
    print >>fileh, '\t'.join(lib_header_fields)
    for fields in lib_metadata:
        print >>fileh, '\t'.join(fields)
    fileh.close()
    # write expr file
    fileh = open(args.expr_file, 'w')
    fields = list(g_header_fields)
    fields.extend(lib_ids)
    print >>fileh, '\t'.join(fields)
    for i in xrange(len(g_inds)):
        fields = list(g_metadata[i])
        fields.extend(map(str, submat[i,:]))
        print >>fileh, '\t'.join(fields)
    fileh.close()
    return 0

if __name__ == '__main__': 
    sys.exit(main())