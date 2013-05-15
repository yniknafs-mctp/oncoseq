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

def read_lines(filename):
    lines = []
    for line in open(filename):
        if not line:
            continue
        if line.startswith("#"):
            continue
        line = line.strip()
        if not line:
            continue
        lines.append(line)
    return lines

def read_table(filename, primary_key='tracking_id', subset=None):
    if subset is None:
        subset = set()
    metadata = []
    inds = []
    keys = []
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    id_ind = header_fields.index(primary_key)
    ind = 0
    for line in fileh:
        fields = line.strip().split('\t')
        key = fields[id_ind]        
        if len(subset) > 0 and key in subset:
            metadata.append(fields)
            keys.append(key)
            inds.append(ind)
        ind += 1
    fileh.close()
    return ind, header_fields, metadata, inds, keys

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--libmetadata", dest="libmetadata", action="store_true", default=False)
    parser.add_argument("--genemetadata", dest="genemetadata", action="store_true", default=False)
    parser.add_argument("library_ids")
    parser.add_argument("gene_ids")
    parser.add_argument("output_file")
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
    library_metadata = []
    library_ids = []
    library_inds = []
    library_table_file = os.path.join(GENE_EXPRESSION_DIR, "library_table.txt")
    fileh = open(library_table_file)
    library_header_fields = fileh.next().strip().split('\t')
    library_id_ind = library_header_fields.index('library_id')
    ncols = 0
    for line in fileh:
        fields = line.strip().split('\t')
        library_id = fields[library_id_ind]
        if library_id in library_id_set:
            if args.libmetadata:
                library_metadata.append(fields)
            library_ids.append(library_id)
            library_inds.append(ncols)
        ncols += 1
    fileh.close()
    if not args.libmetadata:
        library_header_fields = ["library_id"]
    if len(library_inds) == 0:
        print "No libraries found"
        return 1
    # find genes
    gene_metadata = []
    gene_inds = []
    gene_table_file = os.path.join(GENE_EXPRESSION_DIR, "gene_table.txt")
    fileh = open(gene_table_file)
    gene_header_fields = fileh.next().strip().split('\t')
    tracking_id_ind = gene_header_fields.index('tracking_id')
    nrows = 0
    for line in fileh:
        fields = line.strip().split('\t')
        tracking_id = fields[tracking_id_ind]        
        if tracking_id in gene_id_set:
            if not args.genemetadata:
                fields = [tracking_id]
            gene_metadata.append(fields)
            gene_inds.append(nrows)
        nrows += 1
    fileh.close()
    if not args.genemetadata:
        gene_header_fields = ["gene_id"]
    if len(gene_inds) == 0:
        print "No genes found"
        return 1
    # read gene expression
    matrix_file = os.path.join(GENE_EXPRESSION_DIR, "gene_expression.mmap")
    mat = np.memmap(matrix_file, dtype='float32', mode='r', 
                    shape=(nrows,ncols))
    # get subset of matrix
    submat = mat[gene_inds,:]
    submat = submat[:,library_inds]
    # write output
    fileh = open(args.output_file, 'w')
    fields = list(gene_header_fields)
    fields.extend(library_ids)
    print >>fileh, '\t'.join(fields)
    # write library metadata
    ngenemeta = 0
    if len(gene_metadata) > 0:
        ngenemeta = len(gene_metadata[0])
    library_metadata_t = zip(*library_metadata)
    for i in xrange(len(library_metadata_t)):
        fields = ['na'] * (ngenemeta - 1)
        fields.append(library_header_fields[i])
        fields.extend(library_metadata_t[i])
        print >>fileh, '\t'.join(fields)
    # write gene expression
    for i in xrange(len(gene_inds)):
        fields = list(gene_metadata[i])
        fields.extend(map(str, submat[i,:]))
        print >>fileh, '\t'.join(fields)
    fileh.close()
    return 0

if __name__ == '__main__': 
    sys.exit(main())