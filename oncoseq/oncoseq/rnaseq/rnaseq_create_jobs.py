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

# in-place XML prettyprint formatter
def indent_xml(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def generate_analysis_xml(seqdb, output_dir):
    for sample in seqdb.samples.itervalues():
        num_lanes = 0
        for library in sample.libraries:
            num_lanes += sum(1 for lane in library.lanes if lane.qc != "FAIL")
        if num_lanes == 0:
            print "SKIPPED SAMPLE", sample.id
            continue
        output_file = os.path.join(output_dir, "%s.xml" % (sample.id))        
        root = etree.Element("analysis")
        # add patient species, gender age
        elem = etree.SubElement(root, "species")
        elem.text = sample.patient.species
        elem = etree.SubElement(root, "gender")
        elem.text = sample.patient.gender
        elem = etree.SubElement(root, "age")        
        if (sample.patient.age is not None) and sample.patient.age:            
            elem.text = str(int(float(sample.patient.age)))
        else:
            elem.text = ""
        # write sample to xml
        sample_elem = sample.to_xml(root)
        for library in sample.libraries:
            num_lanes = sum(1 for lane in library.lanes if lane.qc != "FAIL")
            if num_lanes == 0:
                print "SKIPPED LIBRARY", sample.id, library.id
                continue
            lib_elem = library.to_xml(sample_elem)
            for lane in library.lanes:
                if lane.qc == "FAIL":
                    print "SKIPPED LANE", lane.id
                    continue
                lane.to_xml(lib_elem)
        f = open(output_file, "w")
        # indent for pretty printing
        indent_xml(root)        
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
    