'''
Created on Jan 4, 2012

@author: mkiyer
'''
import sys

chrominfo = sys.argv[1]
chromsizes = {}
for line in open(chrominfo):
    fields = line.strip().split()
    chrom = fields[0]
    length = int(fields[1])
    chromsizes[chrom] = length

rmsk = sys.argv[2]
for line in open(rmsk):
    if line.startswith("#"):
        continue
    fields = line.strip().split()
    chrom = fields[5]
    start = int(fields[6])
    end = int(fields[7])
    strand = fields[9]
    repname = fields[10]
    repclass = fields[11]
    repfamily = fields[12]
    if repfamily != "rRNA":
        continue
    name = "%s:%d-%d_%s_%s" % (chrom, start, end, strand, repname)
    if chrom == "chrM":
        ensembl_chrom = "MT"
    else:
        ensembl_chrom = chrom[3:]
    ensembl_start = start+1
    ensembl_end = end+1
    if ensembl_chrom not in chromsizes:
        print >>sys.stderr, "could not find chrom", ensembl_chrom
        continue
    print '\t'.join(map(str, [ensembl_chrom, ensembl_start, ensembl_end, strand, name]))