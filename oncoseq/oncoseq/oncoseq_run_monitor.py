'''
Created on Aug 17, 2011

@author: mkiyer
'''
import logging
import argparse
import sqlite3
import os

from oncoseq.lib import rundb
from oncoseq.rnaseq.rnaseq_copy_remote import copy_remote

def monitor_jobs(db, config_file, server_name, output_dir, ssh_port, overwrite):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()
    job_dict = rundb.selectall(c)
    c.close()
    for job_name,status in job_dict.iteritems():
        if status == rundb.STATUS_DONE:
            continue
        logging.info("Monitoring analysis %s current status=%s" % (job_name,status))
        success = copy_remote(output_dir, config_file, server_name, job_name,
                              ssh_port=ssh_port, 
                              overwrite=overwrite)
        if success:
            rundb.setstatus(db, job_name, rundb.STATUS_DONE)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--ssh-port", dest="ssh_port", type=int, default=22)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    parser.add_argument("db_file")
    parser.add_argument("config_file")
    parser.add_argument("server_name")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    # create db
    if not os.path.exists(args.db_file):
        rundb.create_db(args.db_file)
    monitor_jobs(args.db_file,
                 args.config_file,
                 args.server_name,
                 args.output_dir,
                 args.ssh_port,
                 args.overwrite)

if __name__ == '__main__':
    main()
