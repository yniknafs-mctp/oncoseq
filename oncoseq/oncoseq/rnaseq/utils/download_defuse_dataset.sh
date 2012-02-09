wget ftp://ftp.ensembl.org/pub/release-62/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.62.dna.chromosome.*.fa.gz
wget ftp://ftp.ensembl.org/pub/release-62/gtf/homo_sapiens/Homo_sapiens.GRCh37.62.gtf.gz 
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/est.fa.gz 
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/intronEst.txt.gz 
wget ftp://ftp.ncbi.nih.gov/repository/UniGene/Homo_sapiens/Hs.seq.uniq.gz 

# ONE ADDITIONAL STEP MUST BE DONE MANUALLY:
# Navigate to the following URL 
# http://genome.ucsc.edu/cgi-bin/hgTables?command=start. 
# First click at the bottom of the #page where it says To reset all user 
# cart settings (including custom tracks), click here. Then set assembly 
# to Feb. 2009 (GRCh37/hg19), set group to Variation and Repeats, and set 
# track to RepeatMasker. Select gzip compressed and enter an output file 
# at the bottom of the form. Click get output. When the file has finished 
# downloading, gunzip it and set the repeats_filename entry in config.txt 
# to the name of the file.

