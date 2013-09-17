'''
Created on May 22, 2013

@author: mkiyer
'''
import collections
import logging
import argparse
import os
import sys

import xlrd

from oncoseq.rnaseq.lib.libtable import LIBRARY_ROOT_TAG, read_wksheet

def check_library_table_xls(filename):
    """
    checks XLS/XLSX file parameters
    """
    if not os.path.isfile(filename):
        raise OSError("File %s not found or not a regular file" % (filename))
    wkbook = xlrd.open_workbook(filename)
    # check that required sheet names exist
    sheet_names = wkbook.sheet_names()
    if not LIBRARY_ROOT_TAG in sheet_names:
        raise Exception("XLS file missing 'libraries' Sheet")
    if not "parameters" in sheet_names:
        raise Exception("XLS file missing 'parameters' Sheet")
    # read parameters
    params = {"patient": collections.defaultdict(lambda: {}),
              "sample": collections.defaultdict(lambda: {}),
              "library": collections.defaultdict(lambda: {})}
    for field_dict in read_wksheet(wkbook.sheet_by_name("parameters")):
        param_type = field_dict["parameter_type"]
        param_id = field_dict["parameter_id"]
        k = field_dict["parameter_name"]
        v = field_dict["parameter_value"]
        params[param_type][param_id][k] = v
    # read libraries
    patient_ids = set()
    sample_ids = set()
    library_ids = set()
    valid = True
    for field_dict in read_wksheet(wkbook.sheet_by_name(LIBRARY_ROOT_TAG)):
        patient_ids.add(field_dict["patient_id"])
        sample_ids.add(field_dict["sample_id"])
        library_id = field_dict["library_id"]
        if library_id in library_ids:
            logging.error("Found duplicate library id %s" % (library_id))
            valid = False
        library_ids.add(library_id)
    # check parameters
    for k,v in params["patient"].iteritems():
        if k not in patient_ids:
            logging.error("'patient' parameter %s=%s is orphan" % (k,v))
            valid = False
    for k,v in params["sample"].iteritems():
        if k not in sample_ids:
            logging.error("'sample' parameter %s=%s is orphan" % (k,v))
            valid = False
    for k,v in params["library"].iteritems():
        if k not in library_ids:
            logging.error("'library' parameter %s=%s is orphan" % (k,v))
            valid = False
    return valid

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("library_xls_file")
    args = parser.parse_args()
    valid = check_library_table_xls(args.library_xls_file)
    if valid:
        return 0
    return 1

if __name__ == '__main__':
    sys.exit(main())
