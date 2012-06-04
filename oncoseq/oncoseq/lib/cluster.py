'''
Created on Aug 8, 2011

@author: mkiyer
'''
import subprocess
import os
import logging

from config import JOB_SUCCESS, JOB_ERROR

def qstat_user_job_count(remote_address, username, port):
    command = "qstat -u %s | wc -l" % (username)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs: %s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    num_user_jobs = int(res.strip())
    return num_user_jobs

def submit_job_pbs(job_name, 
                   args,
                   num_processors=1,
                   node_processors=1,
                   node_memory=4096,
                   total_memory=None,
                   pbs_script_lines=None, 
                   working_dir=None, 
                   walltime=None, 
                   pmem=None, 
                   mem=None,
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
    pmem: amount of memory allocated to each processor of this job in MB
    mem: amount of total memory to split across all processors (overrides pmem setting)
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
    # ensure number of processors can fit on node
    num_processors = min(num_processors, node_processors)
    resource_fields = ["nodes=%d:ppn=%d" % (1, num_processors)]    
    # setup memory resources
    if mem is not None:
        mem = min(mem, node_memory)
        resource_fields.append("mem=%dmb" % (mem))
    elif pmem is not None:
        max_pmem = float(node_memory) / num_processors
        pmem = min(pmem, max_pmem)
        resource_fields.append("pmem=%dmb" % (pmem))
    else:
        pmem = int(round(float(node_memory) / node_processors, 0))
        mem = pmem * num_processors
        resource_fields.append("mem=%dmb" % (mem))
    # set job walltime
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
                     mem=None,
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



