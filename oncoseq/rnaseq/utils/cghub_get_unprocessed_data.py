'''
Created on Apr 23, 2013

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import xml.etree.cElementTree as etree

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("cghub_xml_file")
    parser.add_argument("processed_aliquots_file")
    args = parser.parse_args()
    # check params
    if not os.path.exists(args.cghub_xml_file):
        parser.error("cghub_xml_file %s not found" % (args.cghub_xml_file))
    if not os.path.exists(args.processed_aliquots_file):
        parser.error("processed_aliquots_file %s not found" % (args.processed_aliquots_file))
    processed_aliquots = set([x.strip() for x in open(args.processed_aliquots_file)])
    # read cghub xml
    tree = etree.parse(args.cghub_xml_file)  
    root = tree.getroot()
    total_results = 0
    unprocessed_aliquots = {}
    for elem in root.findall("Result"):
        total_results += 1
        analysis_id = elem.findtext("analysis_id")
        aliquot_id = elem.findtext("aliquot_id")
        if aliquot_id in processed_aliquots:
            continue
        unprocessed_aliquots[aliquot_id] = analysis_id
    for analysis_id in sorted(unprocessed_aliquots.values()):
        print analysis_id

if __name__ == '__main__':
    sys.exit(main())
