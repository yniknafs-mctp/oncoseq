'''
Created on Nov 2, 2011

@author: mkiyer
'''

def get_read_length(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Sequence length"):
            return int(line.split()[-1])

def get_total_reads(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Total Sequences"):
            return int(line.split()[-1])
