'''
Created on May 8, 2012
@author: mkiyer
'''
import subprocess
import os
import logging
import glob
import datetime
import shutil

from seq import detect_format
from config import JOB_SUCCESS, JOB_ERROR

def scp(src, dst, port, use_compression=False):
    if use_compression:
        args = map(str, ["scp", "-C", "-P", port, src, dst])
    else:
        args = map(str, ["scp", "-P", port, src, dst])
    logging.debug("\targs: %s" % (args))
    return subprocess.call(args)

def ssh_exec(remote_ip, command, port):
    args = map(str, ["ssh", "-p", port, remote_ip, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    return subprocess.call(args) 

def get_remote_file_size(remote_file, remote_address, username, port):
    command = "du -b %s" % (remote_file)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    print args
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    file_size_bytes = int(res.strip().split()[0])
    return file_size_bytes

# TODO: delete this -- redundant implementations (see below) 
#def test_file_exists(remote_address, filename, port):
#    # check that analysis is finished by noting presence of
#    # job completion file
#    command = 'if [ -f "%s" ]; then echo "1"; else echo "0"; fi' % (filename)
#    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
#    logging.debug("\targs=%s" % (args))
#    p = subprocess.Popen(args, stdout=subprocess.PIPE)
#    res = p.communicate()[0]
#    return bool(int(res.strip()))

def test_file_exists(remote_file, remote_address, username, port):
    command = 'test -s %s && echo 1' % (remote_file)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    exists = (res.strip() == "1")
    return exists

def remote_walk(remote_dir, remote_address, username, port):
    command = 'find -P %s -printf "%%P,%%y,%%s\n"' % (remote_dir)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    fileinfo = []
    for line in res.split("\n"):
        fields = line.split(",")
        if (not fields) or (len(fields) <= 1):
            continue
        filename, filetype, filesize = fields
        fileinfo.append((filename, filetype, filesize))
    p.wait()
    return fileinfo

def globus_copy_file(src, dst, remote_address, username, port):
    #
    # I don't know what path format src and dst are in,
    # but they need to be full paths for this
    #
    # Since this isn't technically a remote address (goofy hashtag endpoint) I hardcoded it in
    #
    # As a placeholder until terry adds globus to exds, I put in my credentials on my box -
    #  let me know if you want to test it
    #
    # the port is default (22 i think)
    # 
    #
    logging.debug("Copying files")
    user = "user"
    host = "cli.globusonline.org"
    tx_src  = "user#endpointname"
    tx_dest = "umich#nyx"
    
    args = ["ssh",user + "@" + host,"scp",tx_src + ":" + src,tx_dest + ":" + src];
    retcode = subprocess.call(args)        
    if retcode != 0:
        return JOB_ERROR  
    
    return 0

def copy_to_remote(src, dst, remote_address, username, port, maxsize=(8<<30), tmp_dir="/tmp",login_remote_address=None,):
    # maxsize should be bigger than 1mb
    maxsize = max(maxsize, (1<<20))
    # get compression status of src file
    is_compressed = (detect_format(src) != "txt")
    # check file size and split file if necessary
    src_file_size = os.path.getsize(src)
    if src_file_size > maxsize:
        maxsize_mb = max(1, (maxsize >> 20))
        num_chunks = 1 + ((src_file_size >> 20) / maxsize_mb)
        src_prefix = os.path.basename(os.path.splitext(src)[0])
        # make temp dir for split files
        timestamp_string = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S%f")
        split_dir = os.path.join(tmp_dir, "split_%s" % (timestamp_string))
        os.makedirs(split_dir)
        split_prefix = os.path.join(split_dir, src_prefix)
        # perform file splitting
        logging.debug("Splitting files into %d chunks of size %dmb with file prefix %s" % 
                      (num_chunks, maxsize_mb, split_prefix))
        args = ["split", "-b%dm" % (maxsize_mb), src, split_prefix]
        retcode = subprocess.call(args)        
        if retcode != 0:
            shutil.rmtree(split_dir)
            return JOB_ERROR        
        # copy files to remote location
        logging.debug("Copying files")
        dst_dir = os.path.dirname(dst)
        for split_src_file in glob.glob("%s*" % (split_prefix)):
            split_dst_file = os.path.join(dst_dir, os.path.basename(split_src_file))
            scp(split_src_file, remote_address + ":" + split_dst_file, port, use_compression=(not is_compressed))
        # reconstitute original file from split
        logging.debug("Concatenating chunks")
        command = "cd %s; cat %s* > %s; rm %s*" % (dst_dir, src_prefix, dst, src_prefix)
        retcode = ssh_exec(login_remote_address, command, port)
        # remove local split files
        shutil.rmtree(split_dir)
    else:
        scp(src, remote_address + ":" + dst, port, use_compression=(not is_compressed))
    # check that copy worked by comparing file sizes
    remote_file_size = get_remote_file_size(dst, login_remote_address, username, port)
    if src_file_size != remote_file_size:
        logging.error("\tcopy failed (remote file size != local file size)")
        # cleanup by removing destination file
        retcode = ssh_exec(login_remote_address, "rm %s" % (dst), port)        
        return JOB_ERROR
    else:
        logging.debug("\tremote file size of %d bytes matches local file size" % (src_file_size))
    return 0

def copy_from_remote(src, dst, remote_address, username, port, maxsize=(8<<30), tmp_dir="/tmp"):
    # maxsize should be bigger than 1mb
    maxsize = max(maxsize, (1<<20))
    # get compression status of src file
    is_compressed = (detect_format(src) != "txt")
    # check remote file size and split file if necessary
    src_file_size = get_remote_file_size(src, remote_address, username, port)  
    if src_file_size > maxsize:
        # compute number of chunks to create using split
        maxsize_mb = max(1, (maxsize >> 20))
        num_chunks = 1 + ((src_file_size >> 20) / maxsize_mb)
        src_prefix = os.path.basename(os.path.splitext(src)[0])
        # make remote temp dir for split files
        timestamp_string = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S%f")
        split_dir = os.path.join(tmp_dir, "split_%s" % (timestamp_string))
        logging.debug("Creating tmp dir %s to hold split files" % (split_dir))
        retcode = ssh_exec(remote_address, "mkdir %s" % (split_dir), port)        
        if retcode != 0:
            return JOB_ERROR
        split_prefix = os.path.join(split_dir, src_prefix)
        # perform file splitting
        logging.debug("Splitting files into %d chunks of size %dmb with file prefix %s" % 
                      (num_chunks, maxsize_mb, split_prefix))
        retcode = ssh_exec(remote_address, "split -b%dm %s %s" % (maxsize_mb, src, split_prefix) , port)
        if retcode != 0:
            ssh_exec(remote_address, "rm -rf %s" % (split_dir), port)
            return JOB_ERROR        
        # copy files from remote location
        logging.debug("Copying files")
        dst_dir = os.path.dirname(dst)
        args = ["scp"]
        if not is_compressed:
            args.append("-C")
        args.extend(["-P", port, "-r", remote_address + ":" + split_dir, dst_dir])
        dst_split_dir = os.path.join(dst_dir, os.path.basename(split_dir))
        dst_split_prefix = os.path.join(dst_split_dir, src_prefix)
        retcode = subprocess.call(map(str, args))
        if retcode != 0:
            shutil.rmtree(dst_split_dir)
        # reconstitute original file from split
        logging.debug("Reconstituting original file from split chunks")        
        dstfh = open(dst, "wb")
        for split_src_file in glob.glob("%s*" % dst_split_prefix):
            shutil.copyfileobj(open(split_src_file, "rb"), dstfh)
        dstfh.close()
        # remove local split files
        shutil.rmtree(dst_split_dir)
    else:
        scp(remote_address + ":" + src, dst, port, use_compression=(not is_compressed))
    # check that copy worked by comparing file sizes
    dst_file_size = os.path.getsize(dst)
    if src_file_size != dst_file_size:
        logging.error("\tcopy failed (remote file size != local file size)")
        # cleanup by removing destination file
        os.remove(dst)
        return JOB_ERROR
    else:
        logging.debug("\tremote file size of %d bytes matches local file size" % (src_file_size))
    return JOB_SUCCESS

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    remote_walk('/scratch/arul_flux/med-mctp/projects/rnaseq/runs/job_2012-04-19-150419886085/wcmc_stid581_d', 
                "flux-login.engin.umich.edu", "mkiyer", 22)
    #copy_from_remote("/scratch/arul_flux/med-mctp/projects/rnaseq/runs/lung/lungC/job_2012-04-23-233056054968/lungA85/si_4387/varscan.snvs.txt", 
    #                 "/home/mkiyer/bubba.txt", "flux-login.engin.umich.edu", "mkiyer", 22, maxsize=(1<<20),
    #                 tmp_dir="/scratch/arul_flux/med-mctp/projects/rnaseq/tmp")
