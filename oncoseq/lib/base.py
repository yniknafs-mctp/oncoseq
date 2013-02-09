'''
Created on Aug 3, 2011

@author: mkiyer
'''
import subprocess
import os

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        p = subprocess.Popen([filename], stdout=devnullfh, stderr=devnullfh)
        p.kill()
    except OSError:
        return False
    devnullfh.close()
    return True

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def file_exists_and_nz_size(filename):
    if filename is None:
        return False
    if not os.path.exists(filename):
        return False
    if os.path.getsize(filename) == 0:
        return False
    return True

def parse_bool(s):    
    return True if s[0].lower() == "t" else False

def parse_string_none(s):
    return None if s == "None" else s

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