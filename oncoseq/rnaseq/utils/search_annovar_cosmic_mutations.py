'''
Created on Feb 5, 2013

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import glob

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("mutation_file")
    parser.add_argument("root_dirs", nargs="+")
    args = parser.parse_args()
    # read mutations    
    fileh = open(args.mutation_file)
    header = fileh.next()
    mutdict = {}
    for line in fileh:
        fields = line.strip().split('\t')
        mutdict[fields[0]] = fields
    fileh.close()
    for k,v in mutdict.iteritems():
        logging.debug('Mutation %s %s' % (k,v))
    # read annovar files
    for root_dir in args.root_dirs:
        logging.debug('Root dir: %s' % (root_dir))
        for library_dir in glob.glob(os.path.join(root_dir, "*")):
            if not os.path.isdir(library_dir):
                continue
            logging.debug('Library: %s' % (library_dir))
            library_id = os.path.basename(library_dir)
            annovar_cosmic_file = os.path.join(library_dir, "annovar.hg19_cosmic61_dropped")
            if not os.path.exists(annovar_cosmic_file):
                logging.error('File not found "%s"' % (annovar_cosmic_file))
                continue
            fileh = open(annovar_cosmic_file)
            for line in fileh:
                fields = line.strip().split('\t')
                param_fields = fields[1].split(';')[0]
                valstr = param_fields.split('=')[1]
                #logging.debug('cosmic ids: %s' % (valstr))
                cosmic_ids = valstr.split(',')
                for cosmic_id in cosmic_ids: 
                    if cosmic_id in mutdict:
                        output_fields = [library_id] + mutdict[cosmic_id] + fields 
                        print '\t'.join(output_fields)
            fileh.close()
    

if __name__ == '__main__':
    sys.exit(main())
