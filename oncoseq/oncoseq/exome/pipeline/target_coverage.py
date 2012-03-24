'''
Created on Mar 3, 2012

@author: oabalbin
'''
import logging
import argparse
import subprocess
import sys
import os
import numpy as np
from oncoseq.lib import config

def compute_exon_coverage(bam_file,exon_bed_file,
                   coverage_file, base_qual):
    '''
    It requires samtools and bedtools
    and an bedfile of exon positions
    '''
    args=['samtools','view','-u','-q',str(base_qual),'-f','0x2',bam_file]
    sam_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    args = ['coverageBed','-abam','stdin','-b',exon_bed_file,'-hist']
    f = open(coverage_file, "w")
    retcode = subprocess.call(args, stdin=sam_p.stdout, stdout=f)
    f.close()
    
    if retcode != 0:
        logging.error("samtools pileup failed")
        sam_p.kill()
        if os.path.exists(coverage_file):
            os.remove(coverage_file)
            return config.JOB_ERROR
    return config.JOB_SUCCESS


def summarize_coverage(file,ofile):
    '''
    reads a histogram coverage bed. Summarizes the information
    as chr start end name score strand min max median mean # bases < 10 bp
    chr12    25358179    25362845    KRAS|E0|NM_004985    0    -    32    1    4666    0.0002143
    all     487     3       10718   0.0002799
    '''
    ifile = open(file)
    #of=open(ofile+".add",'w')
    ofcnv=open(ofile,"w")
    poor_th=11
    poor_cov=0
    First=True
    #hd1=["#chr", "start", "end", "name", "score", "strand","min_cov","max_cov","median_cov","mean_cov","bases_with_poor_cov","region_length"]
    hd2=["#name", "chr", "start", "end", "length", "min_coverage", "total_coverage","mean_coverage","bases_with_poor_cov"]
    #of.write(",".join(map(str,hd1)).replace(',','\t')+'\n')
    ofcnv.write(",".join(map(str,hd2)).replace(',','\t')+'\n')
    
    
    for l in ifile:
        
        if l.startswith('all'):
            continue
        
        f=l.strip('\n').split('\t')
        
        name,coverage,bases=f[3],int(f[6]),int(f[7])
            
        if First:
            chr, start, end, name, score, strand, coverage, bases, length,fraction = \
            f[0],f[1],f[2],f[3],f[4],f[5],int(f[6]),int(f[7]),int(f[8]),f[9]
            exon_cov, iname=[coverage*bases],name
            if coverage < poor_th:
                poor_cov += bases
            First=False
            continue
            
        if iname!=name:
            # Print previous record
            total_cov=sum(np.array(exon_cov))
            min_cov,max_cov,median_cov,mean_cov = \
            np.min(exon_cov),np.max(exon_cov),np.median(exon_cov),total_cov/float(length) #np.mean(exon_cov)
            
            #ol=[chr, start, end, name, score, strand,min_cov,max_cov,median_cov,mean_cov,poor_cov,length]
            cnv=[name, chr, start, end, length, min_cov,total_cov,mean_cov, poor_cov]
            
            #of.write(",".join(map(str,ol)).replace(',','\t')+'\n')
            ofcnv.write(",".join(map(str,cnv)).replace(',','\t')+'\n')
            #print iname, name
            exon_cov,poor_cov,iname=[],0,name
            #print exon_cov,poor_cov,iname
            
            
            chr, start, end, name, score, strand, coverage, bases, length,fraction = \
            f[0],f[1],f[2],f[3],f[4],f[5],int(f[6]),int(f[7]),int(f[8]),f[9]
            
            exon_cov,iname=[coverage*bases],name
            if coverage < poor_th:
                poor_cov += bases
            
        else:
            
            exon_cov.append(coverage*bases)
            if coverage < poor_th:
                poor_cov += bases
            
    ifile.close()
    return config.JOB_SUCCESS
    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file")
    parser.add_argument("exon_bed_file")
    parser.add_argument("coverage_file")
    parser.add_argument("coverage_summary_file")
    parser.add_argument("base_quality")
    
    args = parser.parse_args()
    retcode = compute_exon_coverage(args.bam_file, args.exon_bed_file,
                                args.coverage_file,args.base_quality)
    if not retcode:
        return summarize_coverage(args.coverage_file,args.coverage_summary_file)
        
    

if __name__ == '__main__':     
    sys.exit(main())

