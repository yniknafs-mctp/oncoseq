'''
Created on Apr 28, 2012

@author: alebalbin
'''
import os
import sys
import argparse
import logging

import random
import subprocess

def depth_coverage_dilutions(N):
    '''
    return proportions
    '''
    depth_range=range(10,110,10)
    depth_cov_list = [int(N*(i/float(100))) for i in depth_range]
    return depth_cov_list

def fixed_depth_coverage_dilutions():
    '''
    return proportions
    '''
    depth_cov_list=range(10000000,100000000,10000000)
    depth_cov_list.add(5000000)
    return depth_cov_list


def sample_coverage_dilution(fqa,fqb,outputdir,read_basename,N):
    """ 
    get N random headers from a fastq file without reading the
    whole thing into memory 
    python get_subset.py /path/to/reads_1.fastq /path/to/reads_2.fastq 100
    """
    records = sum(1 for _ in open(fqa)) / 4
    #rand_records = sorted([random.randint(0, records - 1) for _ in xrange(N)])
    if N > records:
        msg="The number of reads requested is greater than the total number of reads in the library"
        logging.info(msg)
        N=records

    rand_records = sorted(random.sample(xrange(records), int(N)))
    
    #coverage_dilutions = depth_coverage_dilutions(N)
    coverage_dilutions =fixed_depth_coverage_dilutions()
    ofan,ofbn=read_basename+'_mate0',read_basename+'_mate1'
    
    print N, coverage_dilutions
    for depth_cov in coverage_dilutions:
        # Files for mixed samples
        mix = "_depth%d.fq"%(depth_cov)
        msg="Processing dilution %s"%(mix)
        logging.info(msg)
        print len(rand_records),depth_cov
        
        fha, fhb = open(fqa), open(fqb) 
        ofa,ofb=ofan + mix, ofbn + mix
        ofa,ofb=os.path.join(outputdir,ofa), os.path.join(outputdir,ofb)
        suba, subb = open(ofa, "w"), open(ofb, "w")
        
        # Select number of reads 
        if depth_cov > 0 and depth_cov <= N:
            tumor_records = sorted(random.sample(rand_records, depth_cov ))        
            rec_no = 0
            # Select reads from the tumor sample
            for rr in tumor_records:
                while rec_no < rr:
                    rec_no += 1       
                    for i in range(4): fha.readline()
                    for i in range(4): fhb.readline()
                for i in range(4):
                    suba.write(fha.readline())
                    subb.write(fhb.readline())
                rec_no += 1
        suba.close(), subb.close()
        
        args1=['gzip',ofa]
        r = subprocess.call(args1)
        if r != 0:
            sys.exit(0)
        args2=['gzip',ofb]
        s = subprocess.call(args2)
        if s != 0:
            sys.exit(0)
        
        print >>sys.stderr, "wrote to %s, %s" % (suba.name, subb.name)
        fha.close(), fhb.close()



def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("tumor_read_mate0")
    parser.add_argument("tumor_read_mate1")
    parser.add_argument("outdir")
    parser.add_argument("read_basename")
    parser.add_argument("-n", type=int, dest="number_of_reads", default=100)        
    args = parser.parse_args()
    
    return sample_coverage_dilution(args.tumor_read_mate0, args.tumor_read_mate1, 
                        args.outdir, args.read_basename, args.number_of_reads)

if __name__ == "__main__":
    '''
    It generates a mix of tumor/normal samples for dilution experiments
    to analyze the effect of tumor content on different parameters of a sample.
    '''
    sys.exit(main())
