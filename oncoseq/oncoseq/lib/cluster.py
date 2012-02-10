'''
Created on Aug 8, 2011

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

def qstat_user_job_count(remote_address, username, port):
    command = "qstat -u %s | wc -l" % (username)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    num_user_jobs = int(res.strip())
    return num_user_jobs

def get_remote_file_size(remote_file, remote_address, username, port):
    command = "du -b %s" % (remote_file)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    file_size_bytes = int(res.strip().split()[0])
    return file_size_bytes

def test_file_exists(remote_file, remote_address, username, port):
    command = 'test -s %s && echo 1' % (remote_file)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    exists = (res.strip() == "1")
    return exists

def remote_copy_file(src, dst, remote_address, username, port, maxsize=(8<<30), tmp_dir="/tmp"):
    # maxsize should be bigger than 1mb
    maxsize = max(maxsize, (1<<20))
    # get compression status of src file
    is_compressed = (detect_format(src) != "txt")
    # check file size and split file if necessary
    src_file_size = os.path.getsize(src)
    if src_file_size > maxsize:
        maxsize_mb = max(1, (maxsize >> 20))
        num_chunks = 1 + (src_file_size / maxsize_mb)
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
        retcode = ssh_exec(remote_address, command, port)
        # remove local split files
        shutil.rmtree(split_dir)
    else:
        scp(src, remote_address + ":" + dst, port, use_compression=(not is_compressed))
    # check that copy worked by comparing file sizes
    remote_file_size = get_remote_file_size(dst, remote_address, username, port)
    if src_file_size != remote_file_size:
        logging.error("\tcopy failed (remote file size != local file size)")
        # cleanup by removing destination file
        retcode = ssh_exec(remote_address, "rm %s" % (dst), port)        
        return 1
    else:
        logging.debug("\tremote file size of %d bytes matches local file size" % (src_file_size))
    return 0

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

def submit_job_pbs(job_name, 
                   args,
                   num_processors=1,
                   node_processors=1,
                   node_memory=4096,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   walltime=None, 
                   pmem=None, 
                   deps=None,
                   email=None,
                   stdout_filename=None,
                   stderr_filename=None):
    '''
    job_name: string name of job
    args: list of command-line arguments for job
    num_processors: number of processors to submit with (cannot be greater than the number of processors per node)
    node_processors: number of cores available per node
    node_memory: amount of memory per node (MB)
    pbs_script_lines: list of PBS directives to be added to the script
    working_dir: the "working directory" of the job (allows scripts to access files using relative pathnames)
    walltime: the walltime passed to qsub
    pmem: amount of memory allocated to this job in MB
    deps: 'None' if no dependencies, or a python list of job ids
    email: 'None', or string containing codes 'b', 'a', or 'e' describing when to email
    stdout_filename: string filename for storing stdout
    stderr_filename: string filename for storing stderr
    '''    
    if isinstance(args, basestring):
        args = [args] 
    if pbs_script_lines is None:
        pbs_script_lines = []
    if isinstance(deps, basestring):
        deps = [deps]
    if pmem is None:
        pmem = int(round(float(node_memory)/node_processors, 0))
    # ensure number of processors can fit on node
    num_processors = min(num_processors, node_processors)
    resource_fields = ["nodes=%d:ppn=%d" % (1, num_processors),
                       "pmem=%dmb" % (pmem)]
    if walltime is not None:
        resource_fields.append("walltime=%s" % (walltime))
    # make PBS script
    lines = ["#!/bin/sh",
             "#PBS -N %s" % job_name,
             "#PBS -l %s" % (",".join(resource_fields)),
             "#PBS -V"]
    if email is not None:
        lines.append("#PBS -m %s" % (email))
    if stdout_filename is None:
        stdout_filename = "/dev/null"
    lines.append("#PBS -o %s" % (stdout_filename))
    if stderr_filename is None:
        stderr_filename = "/dev/null"        
    lines.append("#PBS -e %s" % (stderr_filename))
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    lines.extend(pbs_script_lines)
    if working_dir is not None: 
        lines.append("cd %s" % (working_dir))
    lines.append(' '.join(map(str, args)))
    # print to screen
    print '\n'.join(lines)
    # execute script
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p.stdin.write('\n'.join(lines))
    job_id = p.communicate()[0]
    return job_id.strip()

def submit_job_nopbs(job_name, 
                     args,
                     num_processors=1,
                     node_processors=1,
                     node_memory=4096,
                     pbs_script_lines=None, 
                     working_dir=None, 
                     walltime=None, 
                     pmem=None, 
                     deps=None, 
                     email=None,
                     stdout_filename=None,
                     stderr_filename=None):
    curdir = os.getcwd()
    if working_dir is not None:
        os.chdir(working_dir)
    outfh = None
    errfh = None
    if stdout_filename is not None:
        outfh = open(stdout_filename, "w")
    if stderr_filename is not None:
        errfh = open(stderr_filename, "w")
    subprocess.call(' '.join(map(str,args)), shell=True, stdout=outfh, stderr=errfh)
    if stdout_filename is not None:
        outfh.close()
    if stderr_filename is not None:
        errfh.close()
    os.chdir(curdir)



