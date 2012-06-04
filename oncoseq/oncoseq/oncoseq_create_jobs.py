'''
Created on Aug 7, 2011

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import xml.etree.cElementTree as etree

from oncoseq.lib.seqdb import SeqDB
from oncoseq.lib.base import indent_xml

def generate_analysis_xml(seqdb, output_dir):
    for patient in seqdb.patients.itervalues():
        # check for valid data structure
        if not patient.is_valid():
            logging.error("Patient %s not valid" % (patient.id))
            continue
        # add patient        
        root = etree.Element("analysis")
        patient.to_xml(root)
        # indent for pretty printing
        indent_xml(root)
        # write to file
        output_file = os.path.join(output_dir, "%s.xml" % (patient.id))        
        f = open(output_file, "w")
        print >>f, etree.tostring(root)
        f.close()  

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("seqdb_xls_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    seqdb = SeqDB.from_xls(args.seqdb_xls_file)
    generate_analysis_xml(seqdb, args.output_dir)
    
    
if __name__ == '__main__': 
    sys.exit(main())
    