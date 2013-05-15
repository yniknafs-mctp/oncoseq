'''
Created on Apr 30, 2013

@author: mkiyer
'''
import collections

def read_library_table(filename):
    libraries = collections.OrderedDict()
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    library_id_ind = header_fields.index('library_id')
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
        libraries[library_id] = fields
    return header_fields, libraries

def read_gene_metadata_table(filename):
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    tracking_id_ind = header_fields.index('tracking_id')
    genes = collections.OrderedDict()
    for line in fileh:
        fields = line.strip().split('\t')
        tracking_id = fields[tracking_id_ind]        
        genes[tracking_id] = fields
    fileh.close()
    return header_fields, genes