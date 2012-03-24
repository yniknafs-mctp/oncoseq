'''
Created on Dec 8, 2010

@author: mkiyer
'''
import tempfile
import subprocess
import os
import logging
import argparse

def sortIndexSam(input_sam_file, output_bam_file, sort_order='coordinate', picard_path=None, use_mem=2048,):
    '''
    sort sam or bam
    SO=SortOrder Sort order of output file  Required. Possible values: {unsorted, queryname, coordinate} 
    '''
    picard_command = os.path.join(picard_path,"SortSam.jar")

    args = ['java', '-Xmx'+str(use_mem)+'m', '-jar', picard_command, 
            'I=%s' % (input_sam_file), 
            'O=%s' % (output_bam_file),
            'SORT_ORDER=%s' % (sort_order),
            'CREATE_INDEX=true',
            'VALIDATION_STRINGENCY=SILENT'] 
    #'TMP_DIR=/exds/users/oabalbin/tmp/',
    # need to change this TMP_DIR to something else such as a function argument.
    
    retcode = subprocess.call(args)
    
    if retcode != 0:
        raise OSError("command: '%s' returned error code %s" % (' '.join(args), retcode))
    
    return retcode

   
def sam2bam(input_sam_file, output_bam_file, sort_order='coordinate', samtools_path=None, picard_path=None):
    '''
    converts SAM format to BAM format, where reads are sorted
    by position and indexed
    KEY NOTE: Using samtools to convert a sam file to bam allows to control for quality and other parameters
    '''
    
    samtools_bin = 'samtools'
    if samtools_path is not None:
        samtools_bin = os.path.join(samtools_path, samtools_bin)
    try:
        # Define MAPQ to filter for unique and high quality reads
        # It could be passed to the function
        #MAPQ=str(30) # Stringent criteria to select unique reads. 
        MAPQ=str(10)  # Lose criteria to select unique reads, q=1 should be enough.

        # make a temporary file to retain the unsorted BAM
        output_dir = os.path.dirname(output_bam_file)
        fh, tmp_bam_file = tempfile.mkstemp(suffix='bam', prefix='tmp', dir=output_dir)
        os.close(fh)
        #tmp_bam_file=os.path.join(output_dir,'temp.bam')
        #fh = open(tmp_bam_file,'w')
        #fh.close()
        
        args = [samtools_bin, 'view', '-bSq', MAPQ, input_sam_file, '-o', tmp_bam_file]
        print args   
        retcode = subprocess.call(args)
        if retcode != 0:
            raise OSError("command: '%s' returned error code %s" % (' '.join(args), retcode))
        
        else:
            retcode = sortIndexSam(tmp_bam_file, output_bam_file, sort_order, picard_path)
            
    finally:
        if retcode != 0:
            raise OSError("command: sortIndexSam returned error code %s" % (retcode))
        else:
            print "Time to remove the Sam files"
            os.remove(tmp_bam_file)
            os.remove(input_sam_file)
        pass
    
    return retcode
    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--sort-order", dest="sort_order", default="coordinate")
    parser.add_argument("--samtools-path", dest="samtools_path", default=None)
    parser.add_argument("--picard-path", dest="picard_path", default=None)
    parser.add_argument("sam")
    parser.add_argument("bam")
    options = parser.parse_args()
    # Using picard. 
    #return picard_sam2bam(options.sam, options.bam, options.sort_order, options.picard_path)
    return sam2bam(options.sam, options.bam, options.sort_order, options.samtools_path, options.picard_path)

if __name__ == '__main__': 
    import sys
    sys.exit(main())