'''
Created on Aug 16, 2011

@author: mkiyer
'''
import logging
import argparse
import subprocess
import sys
import os
import shutil

from oncoseq.lib import config
from oncoseq.lib.config import AnalysisConfig, PipelineConfig, \
    attach_patient_to_results, validate_patient_results
from oncoseq.lib.remote_server import ssh_exec, scp, test_file_exists, remote_walk, copy_from_remote

# set copy parameters

def copy_remote(output_dir, config_file, server_name, job_name, overwrite):
    # make output directory
    if not os.path.exists(output_dir):
        logging.info("Creating directory: %s" % (output_dir))
        os.makedirs(output_dir)
    # get server parameters        
    pipeline = PipelineConfig.from_xml(config_file)
    server = pipeline.servers[server_name]
    # copy remote analysis file
    logging.info("Copying remote analysis XML file")
    remote_working_dir = os.path.join(server.output_dir, job_name)
    remote_analysis_file = server.address + ":" + os.path.join(remote_working_dir, config.REMOTE_ANALYSIS_XML_FILE)
    analysis_file = os.path.join(output_dir, "analysis_%s.xml" % (job_name)) 
    scp(remote_analysis_file, analysis_file, server.ssh_port) 
    # setup analysis configuration
    analysis = AnalysisConfig.from_xml(analysis_file)
    analysis.attach_to_results(remote_working_dir)
    # remove analysis file
    os.remove(analysis_file)
    # process each sample
    copy_finished = True
    for patient in analysis.patients:
        src_dir = patient.output_dir
        # check for presence of job complete file
        logging.info("Checking analysis for patient %s" % (patient.id))
        if not test_file_exists(server.address, patient.job_complete_file, server.ssh_port):
            logging.error("Patient %s missing 'job.done' status file, job might not be finished" % (patient.id))
            copy_finished = False
            continue
        # prepare destination directory
        dst_dir = os.path.join(output_dir, patient.id)
        logging.info("Checking destination directory: %s" % (dst_dir))
        if os.path.exists(dst_dir):
            if overwrite:
                logging.warning("Overwriting contents of directory: %s" % (dst_dir))
                shutil.rmtree(dst_dir)
            else:
                logging.error("Destination directory exists, cannot copy unless '--overwrite' is set")
                copy_finished = False
                continue
        # walk result directory structure and copy files
        logging.info("Copying results")
        for fileinfo in remote_walk(src_dir, server.address, 
                                    username=None, server.ssh_port):
            filename, filetype, filesize = fileinfo
            src = os.path.join(src_dir, filename)
            dst = os.path.join(dst_dir, filename)
            if filetype == "d":
                # create directory
                logging.debug("Creating directory %s" % (dst))
                os.makedirs(dst)
            else:
                retcode = copy_from_remote(src, dst, server.address, 
                                           username=None, 
                                           port=server.ssh_port, 
                                           maxsize=config.REMOTE_COPY_MAX_SIZE_BYTES,
                                           tmp_dir=server.tmp_dir)
                if retcode != 0:
                    logging.error("Copy of file %s -> %s failed" % (src, dst))
                    copy_finished = False
                    continue
        # validate copy
        logging.info("Validating local copy of results")
        # reattach sample to local results and validate it
        attach_patient_to_results(patient, output_dir)
        is_valid = validate_patient_results(patient)
        if not is_valid:
            logging.error("Copied results are not valid.. deleting")
            copy_finished = False
            shutil.rmtree(dst_dir)
            continue
        # delete remote source directory
        logging.info("Deleting remote data")        
        retcode = ssh_exec(server.address, "rm -rf %s" % (src_dir), server.ssh_port)
        if retcode != 0:
            logging.error("Error deleting sample directory %s" % (src_dir))
            copy_finished = False
    # delete job directory
    if copy_finished:
        logging.info("Deleting remote job")        
        retcode = ssh_exec(server.address, "rm -rf %s" % (remote_working_dir), server.ssh_port)
        if retcode != 0:
            logging.error("Error deleting job directory %s" % (remote_working_dir))
            copy_finished = False
    return copy_finished


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    parser.add_argument("output_dir")
    parser.add_argument("config_file")
    parser.add_argument("server_name")
    parser.add_argument("job_name")
    args = parser.parse_args()
    return copy_remote(args.output_dir,
                       args.config_file,
                       args.server_name,
                       args.job_name,
                       overwrite=args.overwrite)

if __name__ == '__main__': 
    sys.exit(main())