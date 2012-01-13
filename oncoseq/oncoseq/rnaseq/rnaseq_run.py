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

from pantyhose.lib import config
from pantyhose.lib.config import AnalysisConfig, PipelineConfig
from pantyhose.lib.cluster import scp, ssh_exec, qstat_user_job_count
import rnaseq_run_monitor as runmonitor

import pantyhose
_pantyhose_dir = pantyhose.__path__[0]
import pantyhose.pipeline
_pipeline_dir = pantyhose.pipeline.__path__[0] 

def run_local(analysis_file, config_file, server_name, num_processors, 
              db, keep_tmp):
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
    run_cmd = ["python", os.path.join(_pipeline_dir, "run_analysis.py"),
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
               db, keep_tmp):
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
    #
    # make working directory and copy sequence files 
    # 
    for sample in analysis.samples:
        sample_path = os.path.join(server.output_dir, sample.id)
        #
        # make working dir
        #
        logging.info("Making directory for sample %s" % (sample.id))
        ssh_exec(server.address, "mkdir -p %s" % (sample_path), server.ssh_port)
        for library in sample.libraries:
            for lane in library.lanes:
                #
                # copy fastq files to remote location and 
                # replace fastq fields in XML file
                #
                logging.info("Copying lane %s read1 fastq file" % (lane.id))
                ext = os.path.splitext(lane.read1_file)[-1]
                remote_read1_file = os.path.join(sample_path, lane.id + "_1" + ext)
                scp(lane.read1_file, server.address + ":" + remote_read1_file, server.ssh_port)
                lane.read1_file = remote_read1_file
                if lane.read2_file is not None:
                    logging.info("Copying lane %s read2 fastq file" % (lane.id))
                    ext = os.path.splitext(lane.read2_file)[-1]
                    remote_read2_file = os.path.join(sample_path, lane.id + "_2" + ext)
                    scp(lane.read2_file, server.address + ":" + remote_read2_file, server.ssh_port)
                    lane.read2_file = remote_read2_file
    #
    # copy analysis file to remote location
    #
    logging.info("Copying analysis XML file")
    remote_analysis_file = os.path.join(server.output_dir, os.path.basename(analysis_file))         
    analysis.to_xml("tmp_analysis.xml")
    scp("tmp_analysis.xml", server.address + ":" + remote_analysis_file, server.ssh_port)
    os.remove("tmp_analysis.xml")
    #
    # copy configuration file to remote location
    #
    logging.info("Copying pipeline configuration XML file")
    remote_config_file = os.path.join(server.output_dir, os.path.basename(config_file))         
    scp(config_file, server.address + ":" + remote_config_file, server.ssh_port)
    #
    # package and copy source code to remote location
    #    
    logging.info("Packaging source code")
    source_code_file = "tmp_pantyhose_code.tar.gz"
    local_code_py_path = os.path.dirname(pantyhose.__path__[0])
    args = ["tar", "-C", local_code_py_path, "-zcvf", source_code_file, "pantyhose"]
    retcode = subprocess.call(args)
    if retcode != 0:
        logging.error("Error while packaging source code with args: %s" % (" ".join(map(str, args))))
    # create remote directory for code
    timestamp_string = datetime.datetime.now().strftime("%Y-%m-%d-%H%M%S%f")
    remote_code_dir = os.path.join(server.output_dir, "pantyhose_%s" % (timestamp_string)) 
    ssh_exec(server.address, "mkdir -p %s" % (remote_code_dir), server.ssh_port)
    # copy source code
    logging.info("Copying source code")
    remote_code_targz = os.path.join(remote_code_dir, "code.tar.gz")
    scp(source_code_file, server.address + ":" + remote_code_targz, server.ssh_port)
    # unpack source code
    logging.info("Unpacking source code")
    ssh_exec(server.address, "tar -C %s -zxvf %s" % (remote_code_dir, remote_code_targz), server.ssh_port)
    os.remove(source_code_file)
    #
    # submit analysis job
    #
    logging.info("Submitting remote job")
    # setup environment prior to run
    commands = ["source %s" % (server.modules_init_script),
                "cd %s" % (server.output_dir),
                "export PYTHONPATH=%s:$PYTHONPATH" % (remote_code_dir)]                
    for name in pipeline.modules:
        commands.append("module add %s" % (name))
    py_script = os.path.join(remote_code_dir, "pantyhose", "pipeline", "run_analysis.py")
    run_cmd = ["python", py_script, "-p", num_processors, "--rm-fastq"]
    if keep_tmp:
        run_cmd.append("--keep-tmp")
    run_cmd.extend([remote_analysis_file, remote_config_file, server_name])
    commands.append(" ".join(map(str, run_cmd)))
    ssh_exec(server.address, ";".join(commands), server.ssh_port)
    #
    # mark job as "running" in database
    #
    if db is not None:
        if not os.path.exists(db):
            runmonitor.create_db(db)        
        runmonitor.insert(db, analysis_file)
    return 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, dest="num_processors", default=1)
    parser.add_argument("--db", dest="db_file", default=None)
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
                    args.db_file, args.keep_tmp)

if __name__ == '__main__': 
    sys.exit(main())