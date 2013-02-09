'''
Created on Sep 19, 2012

@author: mkiyer
'''
import os
import sys
import argparse
import logging

from oncoseq.lib.seqdb import SeqDB

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("seqdb_xls_file")
    parser.add_argument("lib_table")
    args = parser.parse_args()
    seqdb = SeqDB.from_xls(args.seqdb_xls_file)
    sample_fields = ["sample_type", "disease", "cancer_progression"]
    # attach custom parameters
    sample_params = set()
    for s in seqdb.samples.itervalues():
        sample_params.update(s.params.keys())
    sample_params = sorted(sample_params)
    patient_params = set()
    for p in seqdb.patients.itervalues():
        patient_params.update(p.params.keys())
    patient_params = sorted(patient_params)
    # iterate through library file
    f = open(args.lib_table)
    header_fields = f.next().strip().split('\t')
    header_fields.extend(sample_fields)
    header_fields.extend(sample_params)
    header_fields.extend(patient_params)
    print '\t'.join(header_fields)
    lib_id_col = header_fields.index("library")
    for line in f:
        fields = line.strip().split('\t')
        lib_id = fields[lib_id_col]
        library = seqdb.libraries[lib_id]
        sample = library.sample
        patient = sample.patient
        for field in sample_fields:
            fields.append(getattr(sample, field, None))
        for field in sample_params:
            fields.append(sample.params.get(field))
        for field in patient_params:
            fields.append(patient.params.get(field))
        print '\t'.join(map(str,fields))
                
if __name__ == '__main__': 
    sys.exit(main())
