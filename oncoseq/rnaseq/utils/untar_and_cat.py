'''
Created on Feb 15, 2013

@author: mkiyer
'''
import sys
import os
import subprocess
import shutil
#import glob

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def untar_and_cat(tarfile, tmpdir, output_dir):
    prefix = os.path.basename(tarfile.split('.')[0])
    output_file = os.path.join(output_dir, prefix + ".fastq")    
    print 'prefix', prefix
    print 'output_dir', output_dir
    print 'output_file', output_file
    # untar
    if os.path.exists(tmpdir):
        print 'removing tmpdir', tmpdir
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    args = ['tar', '-zxvf', tarfile, '-C', tmpdir]
    print 'untar args', args
    tarp = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = tarp.communicate()[0]
    files = []
    for line in output.split('\n'):
        line = line.strip()
        filename = os.path.join(tmpdir, line)
        if os.path.exists(filename) and os.path.isfile(filename):
            files.append(os.path.join(tmpdir, line))
    # get files
    files.sort()
    print 'sorted files found', files
    if len(files) == 0:
        return 1
    # concatenate
    cmd = 'cat %s > %s' % (files[0], output_file)
    print 'cat first file args', cmd
    retcode = subprocess.call(cmd, shell=True)
    if retcode == 1:
        if os.path.exists(output_file): os.remove(output_file)
        return 1
    for file in files[1:]:
        cmd = 'cat %s >> %s' % (file, output_file)
        print 'cat append file args', cmd
        retcode = subprocess.call(cmd, shell=True)
        if retcode == 1:
            if os.path.exists(output_file): os.remove(output_file)
            return 1
    # remove tmp dir
    print 'removing tmpdir', tmpdir
    shutil.rmtree(tmpdir)

if __name__ == '__main__':
    if len(sys.argv) < 2:        
        print "no file specified"
        sys.exit(1)
    tarfile = sys.argv[1]
    tmpdir = sys.argv[2]
    outputdir = sys.argv[3]
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    untar_and_cat(tarfile, tmpdir, outputdir)
