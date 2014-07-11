import logging
import argparse
import os
import glob
import sys
from datetime import datetime

def parse_pbs_ppn(pbs_file):
    with open(pbs_file) as f:
        for line in f:
            if line.startswith('#!'):
                continue
            if line.startswith('#PBS'):
                fields = line.strip().split()
                if fields[1] != '-l':
                    continue
                fields = fields[2].split(',')
                for field in fields:
                    if field.startswith('nodes'):
                        nodes,ppn = field.split(':')
                        ppn = int(ppn.split('=')[1])
                        return ppn
    return None

def parse_total_time(stderr_file):
    with open(stderr_file) as f:
        lines = f.readlines()
        timestring = lines[0][1:].split(']')[0]
        t1 = datetime.strptime(timestring, '%a %b %d %H:%M:%S %Z %Y')
        t2 = None
        for line in lines[1:]:
            line = line.strip()
            if line.endswith('Estimating known gene expression with cufflinks'):
                timestring = line[1:].split(']')[0]
                t2 = datetime.strptime(timestring, '%a %b %d %H:%M:%S %Z %Y')
                break
        if t1 is None or t2 is None:
            return None
        return (t2 - t1).total_seconds()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="library_path_file", default=None)
    parser.add_argument("-d", "--dir", dest="input_dirs", action="append", default=None)
    args = parser.parse_args()
    # get library paths either from a file or by searching directories    
    library_paths = []
    if args.library_path_file is not None:
        if not os.path.exists(args.library_path_file):
            parser.error("'%s' not found" % (args.library_path_file))
        for line in open(args.library_path_file):
            path = line.strip()
            if os.path.exists(path) and os.path.isdir(path):
                library_paths.append(path)
            else:
                logging.error("Directory not found at path %s" % (path))
    elif args.input_dirs is not None:
        for input_dir in args.input_dirs:
            if not os.path.exists(input_dir):
                parser.error("input dir '%s' not found" % (input_dir))
            if not os.path.isdir(input_dir):
                parser.error("input dir '%s' not a valid directory" % (input_dir))
            for path in glob.iglob(os.path.join(input_dir, "*")):
                if os.path.exists(path) and os.path.isdir(path):
                    library_paths.append(path)
    else:
        parser.error("Use '-i' to specify a file with a list of result "
                     "paths or '--dir' to search a directory for results")     
    # read libraries
    for output_dir in library_paths:
        logging.debug("Processing %s" % (output_dir))       
        run_pbs_file = os.path.join(output_dir, 'run.pbs')
        stderr_file = os.path.join(output_dir, 'log', 'pbs.stderr')

        if os.path.exists(stderr_file):
            seconds = str(parse_total_time(stderr_file))
        else:
            seconds = 'NA'
        if os.path.exists(run_pbs_file):
            ppn = str(parse_pbs_ppn(run_pbs_file))
        else:
            ppns = 'NA'
        print '\t'.join([output_dir, str(ppn), str(seconds)])

        continue

    
if __name__ == '__main__':
    sys.exit(main())