'''
Created on Jan 17, 2012

@author: mkiyer
'''
import sqlite3

STATUS_RUNNING = "running"
STATUS_DONE = "done"

def create_db(db):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    # Create table
    c.execute('''create table runmonitor (job_name text, status text)''')
    # Save (commit) the changes
    conn.commit()
    c.close()
    
def selectall(c):
    c.execute('select * from runmonitor')
    result_dict = dict((row["job_name"],row["status"]) for row in c)
    return result_dict

def insert(db, job_name):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()    
    job_dict = selectall(c)
    if job_name not in job_dict:
        c.execute('insert into runmonitor values (?,?)', (job_name, STATUS_RUNNING))
    conn.commit()
    c.close()
    
def setstatus(db, job_name, status):
    conn = sqlite3.connect(db)
    conn.row_factory = sqlite3.Row    
    c = conn.cursor()
    job_dict = selectall(c)
    if job_name in job_dict:
        c.execute('update runmonitor set status=? where job_name=?', (status,job_name))
    conn.commit()
    c.close()
