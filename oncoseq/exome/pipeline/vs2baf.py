'''
Created on Feb 6, 2012
@author: oabalbin
Modified on Feb 22, 2012
'''
import sys
import logging
import argparse
import numpy as np
from collections import defaultdict


def vcf_header():
    
    header=['chr\tposition\tcoverage\tbaf']
    
    return header

def read_line(l, fdict,SB_dict):
    '''
    Varscan header for mpileupsomatic calls    
    chrom   position        ref     var     normal_reads1   normal_reads2   normal_var_freq normal_gt       tumor_reads1    tumor_reads2    tumor_var_freq  tumor_gt  
          somatic_status  variant_p_value somatic_p_value tumor_reads1_plus   tumor_reads1_minus       tumor_reads2_plus       tumor_reads2_minus
    '''
    d = defaultdict()
    f = l.strip('\n').replace(':','\t').split('\t')
    #ft= l.strip('\n').split('\t')[SB_dict['SB']]
    for i,v in fdict.iteritems():
        d[i] = f[v]
    #d['SB']=ft
    return d


def print_baf_line(fdict,fields,normal=True):
    '''
    chr    position    coverage    baf
    chr19    323689    68    38
    '''
    # Mandatory fields
    mf = [fields['CHROM'],fields['POS']]
    baf_n=[str(int(fields['DPN1'])+int(fields['DPN2'])),fields['DPN2']]
    baf_t=[str(int(fields['DPT1'])+int(fields['DPT2'])),fields['DPT2']]
    if normal:
        mf=mf+baf_n
    else:
        mf=mf+baf_t
    return mf


def baf_from_varscan(ifile,normal_baf,tumor_baf):
    '''
    Reads the output of varscan snps
    and prints a file in the vcf format
    For exact meaning of the fields see 
    http://varscan.sourceforge.net/using-varscan.html
    '''
    ifile = open(ifile)
    normal_ofile=open(normal_baf,'w')
    tumor_ofile=open(tumor_baf,'w')
    header = ifile.next()
    MAPPING_QUAL=20
    MIN_QUAL=30
    MIN_HET, MAX_HET=0.1, 0.75          
    fdict={'CHROM':0,'POS':1,'REF':2,'ALT':3,
            'DPN1':4,'DPN2':5,'VFN':6,'AAN':7,
            'DPT1':8,'DPT2':9,'VFT':10,'AAT':11,
            'STATUS':12,'VPV':13,'SPV':14}            
    SB_dict={'SBT1P':15,'SBT1M':16,'SBT2P':17,'SBT2M':18}

    normal_ofile.write( ",".join(vcf_header()).replace(',','\n').replace(';',',')+'\n')
    tumor_ofile.write( ",".join(vcf_header()).replace(',','\n').replace(';',',')+'\n')    
    i=0
    for l in ifile:
        #if i < 10:
        fields = read_line(l, fdict,SB_dict)
        # Convert varscan probability to score as SamTools or gatk. Not CAP at 255
        fields['QUAL']= str(np.around(-10 * np.log10( float(fields['SPV'])+1 ),decimals=2)) 
        fields['DP']=str(int(fields['DPN1'])+int(fields['DPN2'])+int(fields['DPT1'])+int(fields['DPT2']))
        #fields['MQ']=str(min(int(fields['MQREF']),fields['MQALT']))
        #fields['SB']=str(round(min(float(fields['SBREF_Plus'])/(float(fields['SBREF_Minus'])+1),
        #                 float(fields['SBALT_Plus'])/(float(fields['SBALT_Minus'])+1)),2) )
        fields['VFN']=str(float(fields['VFN'].replace('%',''))/100.0)
        fields['VFT']=str(float(fields['VFT'].replace('%',''))/100.0)
        fields['FILTER']='PASS' if float(fields['QUAL']) >= MIN_QUAL else 'LowQual'
        fields['ID']='.'#fields['CHROM']+'@'+fields['POS'] 
        
        # Define if the snp is heterozygous or not in the germline.s
        if float(fields['VFN']) >= MIN_HET and float(fields['VFN']) <= MAX_HET:
            #write the normal baf and the tumor baf
            normal_ofile.write( "chr"+",".join(print_baf_line(fdict,fields)).replace(',','\t')+'\n')
            tumor_ofile.write( "chr"+",".join(print_baf_line(fdict,fields,False)).replace(',','\t')+'\n')
        
    normal_ofile.close()
    tumor_ofile.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("snvs_vars_file")
    parser.add_argument("baf_normal_file")
    parser.add_argument("baf_tumor_file")    
    args = parser.parse_args()    

    return baf_from_varscan(args.snvs_vars_file, args.baf_normal_file,args.baf_tumor_file)


if __name__ == '__main__': 
    sys.exit(main())
