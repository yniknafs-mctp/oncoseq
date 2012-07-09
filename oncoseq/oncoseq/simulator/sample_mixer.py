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

def proportion_mixer():
    '''
    return proportions
    '''
    tumor=range(0,110,10)
    normal=range(100,-10,-10)
    return tumor,normal

def sample_mixer(fqa,fqb,bfqa,bfqb,outputdir,read_basename,N):
    """ 
    get N random headers from a fastq file without reading the
    whole thing into memory 
    python get_subset.py /path/to/reads_1.fastq /path/to/reads_2.fastq 100
    """
    records = sum(1 for _ in open(fqa)) / 4
    brecords= sum(1 for _ in open(bfqa)) / 4
    #rand_records = sorted([random.randint(0, records - 1) for _ in xrange(N)])
    if N==None:
        N=min(records,brecords)
        
    rand_records = sorted(random.sample(xrange(records), int(N)))
    brand_records = sorted(random.sample(xrange(brecords), int(N)))
    
    tumor_prop, benign_prop = proportion_mixer()
    ofan,ofbn=read_basename+'_mate0',read_basename+'_mate1'
    
            
    for t, b in zip(tumor_prop, benign_prop):
        # Files for mixed samples
        mix = "_t%d_b%d.fq"%(t,b)
        print mix
        fha, fhb, bha, bhb = open(fqa), open(fqb), open(bfqa), open(bfqb) 
        ofa,ofb=ofan + mix, ofbn + mix
        ofa,ofb=os.path.join(outputdir,ofa), os.path.join(outputdir,ofb)
        suba, subb = open(ofa, "w"), open(ofb, "w")
        
        # Select proportion of tumor/benign 
        if t > 0:
            tumor_records = sorted(random.sample(rand_records, int(len(rand_records)*(t/100.0)) ))        
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
        # Select reads from the benign sample
        if b > 0:
            benign_records = sorted(random.sample(brand_records, int(len(brand_records)*(b/100.0)) ))
            brec_no = 0
            for bb in benign_records:
                while brec_no < bb:
                    brec_no += 1       
                    for i in range(4): bha.readline()
                    for i in range(4): bhb.readline()
                for i in range(4):
                    suba.write(bha.readline())
                    subb.write(bhb.readline())
                brec_no += 1
        suba.close(), subb.close()
        files=[suba,subb]
        args1=['gzip',ofa]
        r = subprocess.call(args1)
        if r != 0:
            sys.exit(0)
        args2=['gzip',ofb]
        s = subprocess.call(args2)
        if s != 0:
            sys.exit(0)
        
        print >>sys.stderr, "wrote to %s, %s" % (suba.name, subb.name)
        fha.close(), fhb.close(), bha.close(), bhb.close()



def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("tumor_read_mate0")
    parser.add_argument("tumor_read_mate1")
    parser.add_argument("normal_read_mate0")
    parser.add_argument("normal_read_mate1")
    parser.add_argument("outdir")
    parser.add_argument("read_basename")
    parser.add_argument("-n", type=int, dest="number_of_reads", default=100)        
    args = parser.parse_args()

    return sample_mixer(args.tumor_read_mate0, args.tumor_read_mate1, args.normal_read_mate0, args.normal_read_mate1, 
                        args.outdir, args.read_basename, args.number_of_reads)

if __name__ == "__main__":
    '''
    It generates a mix of tumor/normal samples for dilution experiments
    to analyze the effect of tumor content on different parameters of a sample.
    '''
    sys.exit(main())
