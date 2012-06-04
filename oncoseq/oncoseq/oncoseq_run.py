'''
Created on Aug 14, 2011

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import datetime
import subprocess

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, PipelineConfig
from oncoseq.lib.cluster import qstat_user_job_count
from oncoseq.lib.remote_server import scp, ssh_exec, test_file_exists, copy_to_remote
from oncoseq.lib.seqdb import FRAGMENT_LAYOUT_PAIRED
from oncoseq.lib import rundb

import oncoseq.pipeline
_pipeline_dir = oncoseq.pipeline.__path__[0] 

def run_local(analysis_file, config_file, server_name, num_processors, 
              db, keep_tmp, local_tmp_dir):
    """
    db: sqlite db file
    """
    #
    # read configuration file
    #
    pipeline = PipelineConfig.from_xml(config_file)
    server = pipeline.servers[server_name]
    #
    # get commands to setup environment prior to run
    #
    commands = ["source %s" % (server.modules_init_script)]
    for name in pipeline.modules:
        basename = name.split("/")[0]
        commands.append("module rm %s" % (basename))
        commands.append("module add %s" % (name))
    #
    # start analysis job
    #
    run_cmd = ["python", os.path.join(_pipeline_dir, "run_oncoseq_analysis.py"),
               "-p", num_processors]
    if keep_tmp:
        run_cmd.append("--keep-tmp")
    run_cmd.extend([analysis_file, config_file, server_name])
    commands.append(" ".join(map(str, run_cmd)))
    command_string = ";".join(commands)
    logging.info("Command string: %s" % (command_string))
    retcode = subprocess.call(command_string, shell=True)
    return retcode

def run_remote(analysis_file, config_file, server_name, num_processors, 
               db, keep_tmp, local_tmp_dir):
    """
    db: sqlite db file
    """
    #
    # read configuration files
    #
    pipeline = PipelineConfig.from_xml(config_file)
    analysis = AnalysisConfig.from_xml(analysis_file)
    # get server parameters
    server = pipeline.servers[server_name]
    # setup remote job directory and job name
    timestamp_string = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S%f")
    job_name = "job_%s" % (timestamp_string)
    job_output_dir = os.path.join(server.output_dir, job_name)
    # update server output directory
    server.output_dir = job_output_dir
    # if running on a PBS queue, check to make sure the job limit is not 
    # exceeded
    if server.pbs:
        #
        # check remote address for user name or try to get
        # from USER environment variable
        #
        if len(server.address.split("@")) > 1:
            username = server.address.split("@")[0]
        else:
            if "USER" in os.environ:
                username = os.environ["USER"]
            else:
                logging.error("unknown username, cannot check server status")
                return config.JOB_ERROR
        #
        # count number of jobs available on cluster
        #
        logging.info("Checking cluster queue status")
        num_user_jobs = qstat_user_job_count(server.address, username, server.ssh_port)
        if num_user_jobs > server.max_user_jobs:
            logging.error("Too many jobs currently running")
            return config.JOB_ERROR
        logging.info("\tjobs in queue: %d" % (num_user_jobs))
    # make working directory for job 
    logging.info("Making directory for job %s" % (job_output_dir))
    ssh_exec(server.address, "mkdir -p %s" % (job_output_dir), server.ssh_port)
    #
    # copy/link to raw sequence files
    #
    for patient in analysis.patients:
        for sample in patient.samples:
            for library in sample.libraries:
                for lane in library.lanes:
                    # handle single versus paired-end sequencing
                    seq_files = ("read1_file",)
                    if lane.fragment_layout == FRAGMENT_LAYOUT_PAIRED:
                        seq_files = ("read1_file", "read2_file")                    
                    for readnum,attrname in enumerate(seq_files):
                        local_file = getattr(lane, attrname)
                        found = False
                        if server.seq_repo_mirror_dir is not None:
                            # test if file exists at remote server sequence 
                            # repository mirror
                            remote_file = os.path.join(server.seq_repo_mirror_dir,
                                                       os.path.basename(local_file))
                            if test_file_exists(remote_file, server.address, 
                                                username, server.ssh_port):
                                logging.info("Found fastq file %s at %s" % (local_file, remote_file))
                                found = True
                                # modify analysis config to point to remote file
                                setattr(lane, attrname, remote_file)
                            else:
                                logging.info("Fastq file %s not found on remote server" % (local_file))
                        if not found:
                            # copy fastq files to remote location and 
                            # replace fastq fields in XML file
                            logging.info("Copying lane %s fastq file" % (lane.id))
                            ext = os.path.splitext(local_file)[-1]
                            remote_file = os.path.join(job_output_dir, "%s_%d%s" % (lane.id, readnum+1, ext))
                            retcode = copy_to_remote(local_file, remote_file, 
                                                     server.address, username, server.ssh_port, 
                                                     maxsize=config.REMOTE_COPY_MAX_SIZE_BYTES,
                                                     tmp_dir=local_tmp_dir)
                            if retcode != 0:
                                logging.error("Copy of read 1 failed")
                                # TODO: cleanup?
                                return config.JOB_ERROR
                            # modify analysis config to point to remote file
                            setattr(lane, attrname, remote_file)
    #
    # copy analysis file to remote location
    #
    logging.info("Copying analysis XML file")
    remote_analysis_file = os.path.join(job_output_dir, config.REMOTE_ANALYSIS_XML_FILE)
    local_analysis_file = os.path.join(local_tmp_dir,"tmp_analysis.xml")
    analysis.to_xml(local_analysis_file)
    scp(local_analysis_file, server.address + ":" + remote_analysis_file, server.ssh_port)
    os.remove(local_analysis_file)
    #
    # copy configuration file to remote location
    #
    logging.info("Copying pipeline configuration XML file")
    remote_config_file = os.path.join(job_output_dir, config.REMOTE_CONFIG_XML_FILE)
    local_config_file = os.path.join(local_tmp_dir, "tmp_pipeline_config.xml")
    pipeline.to_xml(local_config_file)
    scp(local_config_file, server.address + ":" + remote_config_file, server.ssh_port)
    os.remove(local_config_file)    
    #
    # package and copy source code to remote location
    #    
    logging.info("Packaging source code")
    source_code_file = os.path.join(local_tmp_dir, "tmp_oncoseq_code.tar.gz")
    local_code_py_path = os.path.dirname(oncoseq.__path__[0])
    args = ["tar", "-C", local_code_py_path, "-zcvf", source_code_file, "oncoseq"]
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("Error while packaging source code with args: %s" % (" ".join(map(str, args))))
    # copy source code
    logging.info("Copying source code")
    remote_code_targz = os.path.join(job_output_dir, config.REMOTE_CODE_TARGZ_FILE)
    scp(source_code_file, server.address + ":" + remote_code_targz, server.ssh_port)
    # unpack source code
    logging.info("Unpacking source code")
    ssh_exec(server.address, "tar -C %s -zxvf %s" % (job_output_dir, remote_code_targz), server.ssh_port)
    os.remove(source_code_file)
    #
    # submit analysis job
    #
    logging.info("Submitting analysis job")
    # setup environment prior to run
    commands = ["source %s" % (server.modules_init_script),
                "cd %s" % (server.output_dir),
                "export PYTHONPATH=%s:$PYTHONPATH" % (job_output_dir)]                
    for name in pipeline.modules:
        commands.append("module add %s" % (name))
    py_script = os.path.join(job_output_dir, "oncoseq", "pipeline", "run_oncoseq_analysis.py")
    run_cmd = ["python", py_script, "-p", num_processors]
    if keep_tmp:
        run_cmd.append("--keep-tmp")
    run_cmd.extend([remote_analysis_file, remote_config_file, server_name])
    commands.append(" ".join(map(str, run_cmd)))
    ssh_exec(server.address, ";".join(commands), server.ssh_port)
    #
    # mark job as "running" in database
    #
    if db is not None:
        logging.info("Inserting job %s into db %s" % (job_name, db))
        if not os.path.exists(db):
            rundb.create_db(db)
        rundb.insert(db, job_name)
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("Oncoseq version %s" % (oncoseq.__version__))
    logging.info("==========================")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--db", dest="db_file", default=None)
    parser.add_argument("--tmp-dir", dest="local_tmp_dir", help="local temp dir", default="/tmp")
    parser.add_argument("--local", dest="run_local", action="store_true", default=False)
    parser.add_argument("--keep-tmp", action="store_true", dest="keep_tmp", default=False)
    parser.add_argument("analysis_file")
    parser.add_argument("config_file")
    parser.add_argument("server_name")
    args = parser.parse_args()
    if args.run_local:
        run_func = run_local
    else:
        run_func = run_remote    
    return run_func(args.analysis_file, args.config_file, 
                    args.server_name, args.num_processors,
                    args.db_file, args.keep_tmp, args.local_tmp_dir)

if __name__ == '__main__': 
    sys.exit(main())