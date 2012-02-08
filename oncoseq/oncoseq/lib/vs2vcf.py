'''
Created on Feb 6, 2012
@author: oabalbin
'''
import sys
import logging
import argparse
import numpy as np
from collections import defaultdict


def vcf_header():
    '''
    ##fileformat=VCFv4.0
    ##FILTER=<ID=HARD_TO_VALIDATE,Description="MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)">
    ##FILTER=<ID=InDel,Description="Overlaps a user-input mask">
    ##FILTER=<ID=LowQual,Description="Low quality">
    ##FILTER=<ID=STAND_FILTER,Description="QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10">
    ##FILTER=<ID=SnpCluster,Description="SNPs found in clusters">
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
    ##FORMAT=<ID=GL,Number=3,Type=Float,Description="Log-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt; not applicable if site is not biallelic">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AB,Number=1,Type=Float,Description="Allele Balance for hets (ref/(ref+alt))">
    ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
    ##INFO=<ID=Dels,Number=1,Type=Float,Description="Fraction of Reads Containing Spanning Deletions">
    '''
    header=['##fileformat=VCFv4.0','#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE']
    
    return header

def read_line(l, fdict,SB_dict):
    '''
    Chrom   Position        Ref     Var     Cons:Cov:Reads1:Reads2:Freq:P-value     StrandFilter:R1+:R1-:R2+:R2-:pval       SamplesRef      SamplesHet      Sam
    plesHom SamplesNC       Cons:Cov:Reads1:Reads2:Freq:P-value    
    '''
    d = defaultdict()
    f = l.strip('\n').replace(':','\t').split('\t')
    ft= l.strip('\n').split('\t')[SB_dict['SB']]
    for i,v in fdict.iteritems():
        d[i] = f[v]
    d['SB']=ft
    return d

def print_line(fdict,fields):
    '''
    '''
    # Mandatory fields
    mand_f=['ID','CHROM','POS','REF','ALT','QUAL','FILTER'] # ,'FORMAT':'None', 'SAMPLE':'unknown'
    
    mf = [fields['CHROM'],fields['POS'],fields['ID'], 
      fields['REF'], fields['ALT'], fields['QUAL'],
      fields['FILTER']]
    
    opt_f=[]
    for key,val in fields.iteritems():
        if key not in mand_f:
            opt_f.append(key+'='+val)
            
    mf=mf+[",".join(opt_f).replace(',',';')]
    
    return mf


def read_file(ifile,ofile):
    '''
    Reads the output of varscan snps
    and prints a file in the vcf format
    For exact meaning of the fields see 
    http://varscan.sourceforge.net/using-varscan.html
    '''
    ifile = open(ifile)
    ofile=open(ofile,'w')
    header = ifile.next()
    
    MAPPING_QUAL=20
    MIN_QUAL=30
    print header
    '''
    fdict={'CHROM':0,'POS':1,'REF':2,'ALT':3,'DP':5,'DPREF':6,
       'DPALT':7, 'AF':8,'QUAL':9, 'SB':10,'SBREF_Plus':11,
       'SBREF_Minus':12,'SBALT_Plus':13,'SBALT_Minus':14, 'AA':4} 
    '''
    fdict={'CHROM':0,'POS':1,'REF':2,'ALT':3,'DP':5,'DPREF':6,
       'DPALT':7, 'AF':8,'QUAL':9,'AA':4}#, 'SB':10,'SBREF_Plus':11,'SBREF_Minus':12,'SBALT_Plus':13,'SBALT_Minus':14, 'AA':4} 
    SB_dict={'SB':5}

    ofile.write( ",".join(vcf_header()).replace(',','\n')+'\n')
    i=0
    for l in ifile:
        #if i < 10:
        fields = read_line(l, fdict,SB_dict)
        if fields['QUAL']== '0.4999chr6':
            print fields
            continue
        # Convert varscan probability to score as SamTools or gatk. Not CAP at 255
        fields['QUAL']= str(np.around(-10 * np.log10(float(fields['QUAL'])),decimals=2)) 
        #fields['DP']=str(int(fields['DPREF'])+int(fields['DPALT']))
        #fields['MQ']=str(min(int(fields['MQREF']),fields['MQALT']))
        #fields['SB']=str(round(min(float(fields['SBREF_Plus'])/(float(fields['SBREF_Minus'])+1),
        #                 float(fields['SBALT_Plus'])/(float(fields['SBALT_Minus'])+1)),2) )
        fields['AF']=str(float(fields['AF'].replace('%',''))/100.0)
        fields['FILTER']='PASS' if float(fields['QUAL']) >= MIN_QUAL else 'LowQual'
        fields['ID']='.'#fields['CHROM']+'@'+fields['POS']
        fake_sample='\tGT:AD:DP:GL:GQ\t0/0:0,0:0:0,0,0:0'
        #if float(fields['MQREF']) >= MAPPING_QUAL and float(fields['MQALT']) >= MAPPING_QUAL:
        ofile.write( ",".join(print_line(fdict,fields)).replace(',','\t')+fake_sample+'\n')


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("snvs_vars_file")
    parser.add_argument("snvs_vcf_file")
    args = parser.parse_args()    
    return read_file(args.snvs_vars_file, args.snvs_vcf_file)


if __name__ == '__main__': 
    sys.exit(main())
