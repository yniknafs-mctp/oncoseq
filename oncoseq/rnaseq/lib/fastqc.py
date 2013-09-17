'''
Created on Feb 10, 2013

@author: mkiyer
'''
from base import SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT

# FASTQC constants
SANGER_ENCODING = "Sanger / Illumina 1.9"
SOLEXA_ENCODING = "Illumina < 1.3"
ILLUMINA_13_ENCODING = "Illumina 1.3"
ILLUMINA_15_ENCODING = "Illumina 1.5"
ENCODING_VALUES = (SANGER_ENCODING, 
                   SOLEXA_ENCODING, 
                   ILLUMINA_13_ENCODING, 
                   ILLUMINA_15_ENCODING)
ENCODING_TO_QUAL_FORMAT = {SANGER_ENCODING: SANGER_FORMAT,
                           SOLEXA_ENCODING: SOLEXA_FORMAT,
                           ILLUMINA_13_ENCODING: ILLUMINA_FORMAT,
                           ILLUMINA_15_ENCODING: ILLUMINA_FORMAT}

def get_most_common_sequence_length(fastqc_data_file):
    fileh = open(fastqc_data_file)
    for line in fileh:
        if line.startswith(">>Sequence Length Distribution"):
            break
    fileh.next()
    most_common_length = None
    most_common_count = 0
    for line in fileh:
        if line.startswith("#"):
            continue
        if line.startswith(">>END_MODULE"):
            break
        fields = line.strip().split('\t')
        length = fields[0]
        count = float(fields[1])
        if (count >= most_common_count):
            most_common_length = length
            most_common_count = count
    fileh.close()
    lengths = map(int, most_common_length.split('-'))
    return int(round(float(sum(lengths)) / len(lengths)))

def get_sequence_length(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Sequence length"):
            return int(line.split()[-1])

def get_total_sequences(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Total Sequences"):
            return int(line.split()[-1])

def get_encoding(fastqc_data_file):
    for line in open(fastqc_data_file):
        if not line: continue
        line = line.strip()
        if line.startswith("Encoding"):
            return line.split("\t")[-1]