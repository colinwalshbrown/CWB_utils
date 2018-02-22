#!/usr/bin/env python
"""
xgrid_tools.py
"""
import pickle
import tempfile
import os
import sys
import re
from subprocess import Popen,PIPE
from time import sleep

XGRID_HOST = "kermit.QB3.Berkeley.EDU"
ULIMIT = 20

class xgrid_runner(object):
    def __init__(self,jobs=None):
        self.jobs = jobs
        self.completed_jobs = []
        self.failed_jobs = []
        self.all_done = False

    def add_job(self,job):
        if self.jobs:
            self.jobs.append(job)
        else:
            self.jobs = [job]

    def run_jobs(self,hold=1,verbose=False):
        self.tempjoblist = []
        self.overflowjoblist = []
        if len(self.jobs) > ULIMIT:
            self.tempjoblist = self.jobs[0:ULIMIT]
            self.overflowjoblist = self.jobs[ULIMIT:]
            
        else:
            self.tempjoblist = self.jobs
            
        [x.start() for x in self.tempjoblist]
        self.completed_jobs = []
        if (verbose):
            print [(x.jobID,x.obj) for x in self.tempjoblist]
        while(not self.all_done):
            self.donejobs = [x for x in self.tempjoblist if ((x.checkstatus() == "DONE") or (x.checkstatus() == "FAIL"))]
            self.completed_jobs += [x for x in self.tempjoblist if ((x.checkstatus() == "DONE") and (x not in self.completed_jobs))]
            self.failed_jobs += [x for x in self.tempjoblist if (x.checkstatus() == "FAIL")]
            remaining = len(self.tempjoblist) - len(self.donejobs) 
            if verbose:
                print "Failed: %i\tCompleted: %i\tRunning: %i\tIn Queue: %i" % (len(self.failed_jobs),len(self.completed_jobs),remaining,len(self.overflowjoblist))
            if (remaining < ULIMIT) and len(self.overflowjoblist) > 0:
                self.runningjobs = [x for x in self.tempjoblist if ((x.checkstatus() != "DONE") and (x.checkstatus() != "FAIL"))]
                if (len(self.overflowjoblist) > (ULIMIT-remaining)):
                    self.newjobs = self.overflowjoblist[0:(ULIMIT-remaining)]
                    self.overflowjoblist = self.overflowjoblist[(ULIMIT-remaining):]
                else:
                    self.newjobs = self.overflowjoblist
                    self.overflowjoblist = []
                [x.start() for x in self.newjobs]
                self.tempjoblist = self.runningjobs + self.newjobs
            elif (remaining == 0):
                self.all_done = True
                [x.cleanup() for x in self.tempjoblist]
            [x.cleanup() for x in self.donejobs]
            sleep(hold)
        
    def get_done(self):
        if not self.all_done:
            print >> sys.stderr, "Warning: not all jobs have been completed"
        return self.completed_jobs

class xgrid_job(object):
    def __init__(self,run_cmd=None,dir=None,host=None,name=None): 
        self.run_cmd = run_cmd
        self.dir = dir
        if host:
            self.host = host
        else:
            self.host = XGRID_HOST
        self.jobID = None
        self.started = False
        self.done = False
        self.failed = False
        self.name = name
        object.__init__(self)

    def start(self):
        xgrid_args = ["xgrid",
                      "-h",
                      XGRID_HOST,
                      "-auth",
                      "Kerberos",
                      "-job",
                      "submit",
                      "-in",
                      self.dir,
                      self.run_cmd]
        print xgrid_args
        process = Popen(xgrid_args,stdout=PIPE)
        output = process.communicate()[0]
        job_search = re.search("jobIdentifier = (\S+);",output)
        if not job_search:
            print xgrid_args
            print output
            print self.run_cmd
        self.jobID = job_search.groups()[0]
        self.started = True
        return self.jobID

    def checkstatus(self):
        xgrid_status = ['xgrid',
                        '-h',
                        XGRID_HOST,
                        '-auth',
                        'Kerberos',
                        '-job',
                        'attributes',
                        '-id',
                        self.jobID]
        print xgrid_status,
        print " " + self.name
        status_get = Popen(xgrid_status,stdout=PIPE).communicate()[0]
        status = re.search("jobStatus = (\S+);",status_get)
        if status and status.groups()[0] == "Finished":
            if not self.done_obj:
                self.done_obj = self.get_results()
            
            if self.done_obj:
                self.done = True
                return "DONE"
            else:
                self.failed = True
                return "FAIL"
        
        elif status and status.groups()[0] == "Failed":
            self.failed = True
            return "FAIL"
        elif status and status.groups()[0] == "Started":
            return "RUNNING"
        else:
            return None

    def get_results(self):
        xgrid_result = ["xgrid",
                        "-h",
                        XGRID_HOST,
                        "-auth",
                        "Kerberos",
                        "-job",
                        "results",
                        "-id",
                        self.jobID]
        return Popen(xgrid_result,stdout=PIPE).communicate()[0]

    def cleanup(self):
        self.delete()

    def delete(self):
        if self.started:
            delete_args = ["xgrid",
                           "-h",
                           XGRID_HOST,
                           "-auth",
                           "Kerberos",
                           "-job",
                           "delete",
                           "-id",
                           self.jobID]
            Popen(delete_args).communicate()
    
    
    

class xgrid_object_job(object):
    
    def __init__(self,obj,host=None,name=None):
        self.obj = obj
        self.done_obj = None

        if host:
            self.host = host
        else:
            self.host = XGRID_HOST

        self.start_pickle_file = None
        self.done_pickle_file = None
        self.script_file = None
        self.jobID = None
        self.started = False
        self.done = False
        self.failed = False
        self.name = name
        object.__init__(self)

    def start(self,runfn='runXgrid',check=1):
        # create tempfiles, one for the object and one for the script
        self.start_pickle_file = tempfile.NamedTemporaryFile()
        self.start_pickle_name = self.start_pickle_file.name.split("/")[-1]
        self.done_pickle_file = tempfile.NamedTemporaryFile()
        self.done_pickle_name = self.done_pickle_file.name.split("/")[-1]
        self.script_file = tempfile.NamedTemporaryFile(suffix=".py")
        self.tempdir = tempfile.tempdir
        
        # dump the object into pickle_file
        start_pickle = pickle.dumps(self.obj)
        print >> self.start_pickle_file, start_pickle,
        self.start_pickle_file.flush()

        # write the temp script to script_file
        script_text = self.write_script(runfn)
        print >> self.script_file, script_text
        self.script_file.flush()
        print self.script_file.name,
        # make script_file executable
        os.chmod(self.script_file.name,0777)

        # start script_file
        script_file_nodir = self.script_file.name.split("/")[-1]
        xgrid_args = ["xgrid",
                      "-h",
                      XGRID_HOST,
                      "-auth",
                      "Kerberos",
                      "-job",
                      "submit",
                      "-in",
                      self.tempdir,
                      script_file_nodir]
        process = Popen(xgrid_args,stdout=PIPE)
        output = process.communicate()[0]
        job_search = re.search("jobIdentifier = (\S+);",output)
        if not job_search:
            print output
            print self.obj
        self.jobID = job_search.groups()[0]
        self.started = True
        return self.jobID

    def checkstatus(self):
        xgrid_status = ['xgrid',
                        '-h',
                        XGRID_HOST,
                        '-auth',
                        'Kerberos',
                        '-job',
                        'attributes',
                        '-id',
                        self.jobID]
        print xgrid_status,
        print " " + self.name
        status_get = Popen(xgrid_status,stdout=PIPE).communicate()[0]
        #print status_get
        status = re.search("jobStatus = (\S+);",status_get)
        #print status.groups()
        if status and status.groups()[0] == "Finished":
            if not self.done_obj:
                self.done_obj = self.get_results()
            
            if self.done_obj:
                self.done = True
                return "DONE"
            else:
                self.failed = True
                return "FAIL"
        
        elif status and status.groups()[0] == "Failed":
            self.failed = True
            return "FAIL"
        else:
            return None

    def get_results(self):
        xgrid_result = ["xgrid",
                        "-h",
                        XGRID_HOST,
                        "-auth",
                        "Kerberos",
                        "-job",
                        "results",
                        "-id",
                        self.jobID]
        result = Popen(xgrid_result,stdout=PIPE).communicate()[0]
        err = 0
        #xgrid_err = ["xgrid",
        #             "-h",
        #             XGRID_HOST,
        #             "-auth",
        #             "Kerberos",
        #             "-job",
        #             "results",
        #             "-se",
        #             "-id",
        #             self.jobID]
        #err = Popen(xgrid_err,stdout=PIPE).communicate()[0]
        try:
            return pickle.loads(result)
        except:
            print "job %s failed:" % self.jobID
            return None

    def cleanup(self):
        if self.start_pickle_file:
            self.start_pickle_file.close()

        if self.done_pickle_file:
            self.done_pickle_file.close()

        if self.script_file:
            self.script_file.close()

    def delete(self):
        if self.started():
            delete_args = ["xgrid",
                           "-h",
                           XGRID_HOST,
                           "-auth",
                           "Kerberos",
                           "-job",
                           "delete",
                           "-id",
                           self.jobID]
            Popen(delete_args).communicate()
        
    def write_script(self,runfn):
        script_text = \
"""#!/usr/bin/env python\n
\"""
Created by Colin 'Fucking' Brown's xgrid_Runner. 
If you come across this file in error, rest assured it is YOUR FAULT.  Back away slowly, and speak NOT A WORD OF THIS to ANYONE.
\"""
import pickle
import sys        

def _main():
\tstart_file = open('%s')
\tdone_file = open('%s','w+b')
\tsys.path = %s
\ttry:
\t\tobj = pickle.load(start_file)
\t\tobj.%s()
\t\tprint pickle.dumps(obj)
\tfinally: 
\t\tstart_file.close()
\t\tdone_file.close()
if __name__ == "__main__":
\t_main()
""" %  (self.start_pickle_name,self.done_pickle_name,sys.path,runfn)
        return script_text
