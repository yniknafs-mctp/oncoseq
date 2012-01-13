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

from pantyhose.lib import config
from pantyhose.lib.config import AnalysisConfig, attach_sample_to_results, validate_sample_results
from pantyhose.lib.cluster import ssh_exec

def test_file_exists(remote_address, filename, port):
    # check that analysis is finished by noting presence of
    # job completion file
    command = 'if [ -f "%s" ]; then echo "1"; else echo "0"; fi' % (filename)
    args = map(str, ["ssh", "-p", port, remote_address, '%s' % (command)])
    logging.debug("\targs=%s" % (args))
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    return bool(int(res.strip()))

def copy_remote(analysis_file, output_dir,  remote_address, 
                remote_working_dir, ssh_port, overwrite):
    # make output directory
    if not os.path.exists(output_dir):
        logging.info("Creating directory: %s" % (output_dir))
        os.makedirs(output_dir)
    #
    # read configuration files
    #
    analysis = AnalysisConfig.from_xml(analysis_file)
    analysis.attach_to_results(remote_working_dir)
    # process each sample
    copy_finished = True
    for sample in analysis.samples:
        src_dir = sample.output_dir
        # check that analysis is finished by noting presence of
        # job completion file
        logging.info("Checking that analysis is complete")
        if not test_file_exists(remote_address, sample.job_complete_file, ssh_port):
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
        args = map(str, ["scp", "-P", ssh_port, "-r", remote_address + ":" + src_dir, dst_dir])
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
        retcode = ssh_exec(remote_address, "rm -rf %s" % (src_dir), ssh_port)
        if retcode != 0:
            logging.error("Error deleting source directory %s" % (src_dir))
            copy_finished = False
    return copy_finished

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--ssh-port", dest="ssh_port", type=int, default=22)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    parser.add_argument("analysis_file")
    parser.add_argument("output_dir")
    parser.add_argument("remote_address")
    parser.add_argument("remote_working_dir")
    args = parser.parse_args()
    return copy_remote(args.analysis_file, 
                       args.output_dir,
                       args.remote_address,
                       args.remote_working_dir,
                       args.ssh_port,
                       overwrite=args.overwrite)
    
if __name__ == '__main__': 
    sys.exit(main())