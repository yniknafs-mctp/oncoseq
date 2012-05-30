'''
Created on May 27, 2011

@author: oabalbin

This script converts the quality score of the a fq file
from illumina to sanger, so it can be used in general 
for gatk, bwa, and samtools. Note that after casava 1.8
illumina is going to produce sanger quality score by default.
'''

#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord
import sys
from optparse import OptionParser
from Bio import SeqIO
from base import JOB_SUCCESS, JOB_ERROR
import subprocess

def illu2sanger(ifile,ifile_fq,file_report):
    '''
    convert the file with illumina quality
    to fastq-sanger qualitiy
    '''
    count = SeqIO.convert(ifile, "fastq-illumina", ifile_fq, "fastq-sanger")
    if count > 0:
        print file_report
        JOB_SUCCESS_FILE = open(file_report,'w')
        JOB_SUCCESS_FILE.write(str(count))
        JOB_SUCCESS_FILE.close()
        subprocess.Popen(['mv',ifile_fq,ifile])
        return JOB_SUCCESS
    else:
        return JOB_ERROR
    
    
if __name__ == '__main__':
    optionparser = OptionParser("usage: %prog [options] ")
    
    optionparser.add_option("-i", "--ifile", dest="ifile",
                            help="file of reads with quality in a given format")
    optionparser.add_option("-o", "--ifile_fq", dest="ifile_fq",
                            help="file of reads with qualities in sanger scale") 
    optionparser.add_option("-p", "--platform", dest="platform",
                            help="platform that determines scale, solexa, illumina, sanger") 
    optionparser.add_option("-r", "--report_file", dest="report_file",
                            help="report file") 


    (options, args) = optionparser.parse_args()
    if options.platform=="illumina":
        sys.exit(illu2sanger(options.ifile,options.ifile_fq, options.report_file))