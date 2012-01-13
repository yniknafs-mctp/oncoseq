'''
Created on Aug 17, 2011

@author: mkiyer
'''
import logging
import argparse
import sqlite3
import os

from rnaseq_copy_remote import copy_remote

STATUS_RUNNING = "running"
STATUS_DONE = "done"

def create_db(db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    # Create table
    c.execute('''create table runmonitor (analysis_file text, status text)''')
    # Save (commit) the changes
    conn.commit()
    c.close()

def selectall(c):
    c.execute('select * from runmonitor')
    files = dict((row["analysis_file"],row["status"]) for row in c)
    return files

def insert(db, analysis_file):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()    
    files = selectall(c)
    if analysis_file not in files:
        c.execute('insert into runmonitor values (?,?)', (analysis_file, STATUS_RUNNING))
    conn.commit()
    c.close()

def setstatus(db, analysis_file, status):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()    
    files = selectall(c)
    if analysis_file in files:
        c.execute('update runmonitor set status=? where analysis_file=?', (status,analysis_file))
    conn.commit()
    c.close()

def monitor_jobs(db, output_dir, remote_address, remote_working_dir, ssh_port, overwrite):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()
    files = selectall(c)
    c.close()
    for f,status in files.iteritems():
        if status == STATUS_DONE:
            continue
        logging.info("Monitoring analysis %s current status=%s" % (f,status))
        success = copy_remote(f, output_dir,
                              remote_address,
                              remote_working_dir,
                              ssh_port=ssh_port,
                              overwrite=overwrite)
        if success:
            setstatus(db, f, STATUS_DONE)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--ssh-port", dest="ssh_port", type=int, default=22)
    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False)
    parser.add_argument("db_file")
    parser.add_argument("output_dir")
    parser.add_argument("remote_address")
    parser.add_argument("remote_working_dir")
    args = parser.parse_args()
    # create db
    if not os.path.exists(args.db_file):
        create_db(args.db_file)
    monitor_jobs(args.db_file,
                 args.output_dir,
                 args.remote_address,
                 args.remote_working_dir,
                 args.ssh_port,
                 args.overwrite)

if __name__ == '__main__':
    main()
