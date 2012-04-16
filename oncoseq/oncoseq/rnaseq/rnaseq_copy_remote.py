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
from oncoseq.lib.config import AnalysisConfig, PipelineConfig, attach_sample_to_results, validate_sample_results
from oncoseq.lib.cluster import ssh_exec, scp

def test_file_exists(remote_address, filename, port):
    # check that analysis is finished by noting presence of
    # job completion file
    command = 'if [ -f "%s" ]; then echo "1"; else echo "0"; fi' % (filename)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs=%s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    return bool(int(res.strip()))

def copy_remote(output_dir, config_file, server_name, job_name,
                ssh_port, overwrite):
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
    for patient, patient_samples in analysis.patients.iteritems():
        src_dir = patient.output_dir
        for sample in patient_samples:
            #src_dir = sample.output_dir
            # check that analysis is finished by noting presence of
            # job completion file
            logging.info("Checking that analysis is complete")
            if not test_file_exists(server.address, sample.job_complete_file, ssh_port):
                logging.error("Sample %s missing 'job.done' status file, job might not be finished" % (sample.id))
                copy_finished = False
                continue
            # prepare destination directory
            dst_dir = os.path.join(output_dir, sample.id)
            logging.info("Preparing destination directory: %s" % (dst_dir))
            if os.path.exists(dst_dir):
                if overwrite:
                    logging.warning("Overwriting contents of directory: %s" % (dst_dir))
                    shutil.rmtree(dst_dir)
                else:
                    logging.error("Destination directory exists, cannot copy unless '--overwrite' is set")
                    copy_finished = False
                    continue
            # copy to destination directory
            logging.info("Copying results")
            args = map(str, ["scp", "-P", ssh_port, "-r", server.address + ":" + src_dir, dst_dir])
            retcode = subprocess.call(args)
            if retcode != 0:
                logging.error("Error copying results from %s to %s" % (src_dir, dst_dir))
                copy_finished = False
                continue
            # validate copy
            logging.info("Validating copy")
            # reattach sample to local results and validate it
            attach_sample_to_results(sample, output_dir)
            is_valid = validate_sample_results(sample)
            if not is_valid:
                logging.error("Copied sample results are not valid.. deleting")
                copy_finished = False
                shutil.rmtree(dst_dir)
                continue
            # delete remote source directory
            logging.info("Deleting remote data")        
            retcode = ssh_exec(server.address, "rm -rf %s" % (src_dir), ssh_port)
            if retcode != 0:
                logging.error("Error deleting sample directory %s" % (src_dir))
                copy_finished = False
    # delete job directory
    if copy_finished:
        logging.info("Deleting remote job")        
        retcode = ssh_exec(server.address, "rm -rf %s" % (remote_working_dir), ssh_port)
        if retcode != 0:
            logging.error("Error deleting job directory %s" % (remote_working_dir))
            copy_finished = False
    return copy_finished


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--ssh-port", dest="ssh_port", type=int, default=22)
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
                       args.ssh_port,
                       overwrite=args.overwrite)

if __name__ == '__main__': 
    sys.exit(main())