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
    
    header=['##fileformat=VCFv4.0',
    '##FILTER=<ID=LowQual;Description="Low quality Somatic Qual < 30">',
    '##INFO=<ID=DPN1;Number=Depth of coverage supporting the reference allele in the normal sample">',
    '##INFO=<ID=DPN2;Number=Depth of coverage supporting the alternative allele in the normal sample">',
    '##INFO=<ID=VFN;Number=Variant fraction in the normal sample">',
    '##INFO=<ID=AAN;character=Expected AA in the normal sample">',
    '##INFO=<ID=DPT1;Number=Depth of coverage supporting the reference allele in the tumor sample">',
    '##INFO=<ID=DPT2;Number=Depth of coverage supporting the alternative allele in the tumor sample">',
    '##INFO=<ID=VFT;Number=Variant fraction in the tumor sample">',
    '##INFO=<ID=AAT;character=Expected AA in the tumor sample">',
    '##INFO=<ID=STATUS;character=Germline; Somatic; or LOH">',
    '##INFO=<ID=VPV;Number=pvalue for the variant">',
    '##INFO=<ID=VPV;Number=pvalue for the somatic variant">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE']
    
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


def convert_varscan_to_vcf(ifile,ofile):
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
              
    fdict={'CHROM':0,'POS':1,'REF':2,'ALT':3,
            'DPN1':4,'DPN2':5,'VFN':6,'AAN':7,
            'DPT1':8,'DPT2':9,'VFT':10,'AAT':11,
            'STATUS':12,'VPV':13,'SPV':14}            
    SB_dict={'SBT1P':15,'SBT1M':16,'SBT2P':17,'SBT2M':18}

    ofile.write( ",".join(vcf_header()).replace(',','\n').replace(';',',')+'\n')
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
        fake_sample='\tGT:AD:DP:GL:GQ\t0/0:0,0:0:0,0,0:0'
        #if float(fields['MQREF']) >= MAPPING_QUAL and float(fields['MQALT']) >= MAPPING_QUAL:
        ofile.write( ",".join(print_line(fdict,fields)).replace(',','\t')+'\n') #+fake_sample+'\n')


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("snvs_vars_file")
    parser.add_argument("snvs_vcf_file")
    args = parser.parse_args()    
    #snvs_vcf_file = args.snvs_vars_file.replace('.vcf','.2.vcf')

    return convert_varscan_to_vcf(args.snvs_vars_file, args.snvs_vcf_file)


if __name__ == '__main__': 
    sys.exit(main())
