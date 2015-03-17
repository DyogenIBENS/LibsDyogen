#!/usr/bin/python
# -*- coding: utf-8 -*-

# author:  Adrian Sampson (sampsyo) http://homes.cs.washington.edu/~asampson/
# Found on https://github.com/sampsyo/clusterfutures

"""
Abstracts access to a Condor cluster via its command-line tools.
"""

import subprocess
import re
import os
import threading
import time

LOG_FILE = "condorpy.log"
OUTFILE_FMT = "condorpy.stdout.%s.log"
ERRFILE_FMT = "condorpy.stderr.%s.log"

def call(command, stdin=None):
    """Invokes a shell command as a subprocess, optionally with some
    data sent to the standard input. Returns the standard output data,
    the standard error, and the return code.
    """
    if stdin is not None:
        stdin_flag = subprocess.PIPE
    else:
        stdin_flag = None
    proc = subprocess.Popen(command, shell=True, stdin=stdin_flag,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate(stdin)
    return stdout, stderr, proc.returncode

class CommandError(Exception):
    """Raised when a shell command exits abnormally."""
    def __init__(self, command, code, stderr):
        self.command = command
        self.code = code
        self.stderr = stderr

    def __str__(self):
        return "%s exited with status %i: %s" % (repr(self.command),
                                                 self.code,
                                                 repr(self.stderr))

def chcall(command, stdin=None):
    """Like ``call`` but raises an exception when the return code is
    nonzero. Only returns the stdout and stderr data.
    """
    stdout, stderr, code = call(command, stdin)
    if code != 0:
        raise CommandError(command, code, stderr)
    return stdout, stderr

def submit_text(job):
    """Submits a Condor job represented as a job file string. Returns
    the cluster ID of the submitted job.
    """
    out, _ = chcall('condor_submit -v', job)
    jobid = re.search(r'Proc (\d+)\.0', out).group(1)
    return int(jobid)

def submit(script, universe="vanilla", arguments=None, mail=None,  # group name of the job
           jobGroup='<group>', niceUser=False, requirements='(Memory > 1024)', priority=0,  # maximum number of simultaneous jobs with the same group name
           maxSimultaneousJobsInGroup=100, log=LOG_FILE, outfile=OUTFILE_FMT % "$(Cluster)", errfile=ERRFILE_FMT % "$(Cluster)"):

    """Starts a Condor job based on specified parameters. A job
    description is generated. Returns the cluster ID of the new job.

    # examples of options:
    # 1) requirements:
    # requirements = (Memory > 1024) && (Machine != "bioclust01.bioclust.biologie.ens.fr") && (slotid <= 6)
    # (slotid <= 6) means that the job won't occupy more than 6 thread on a multi-core machine.
    # This option may be interesting if the job uses several threads on a same machine.
    # For instance Blast uses 4 threads and it would be interesting to send the job with (slotid <= 4) on an octo-core
    # machine.
    # 2) priority:
    # Each job may have a priority level specified by 'priority', an int >= 0.
    # Condor executes jobs by decreasing priority levels.
    # If a job1 has a priority1 and job2 has a priority2 and if priority1 > priority2, job1 will be executed before job2
    # This won't change the priority between two different users, it just changes the order of jobs of a same user.
    # 3) niceUser:
    # If the user is planning to send a big amount of jobs that do not take long to execute, he may execute them with
    # the niceUser option set to True (and no limit of simultaneous runs), thus all other users will have the priority
    # over him.
    # 4) maxSimRuns is an upper limit of simultaneous runs on all the machines of the cluster. The user can send an
    # unlimited amount of jobs that will be queued but no more than 100 jobs will be executed simultaneously.
    # This limit limit the saturation of cores and keep some cores free for other users.
    # This option is especially important if jobs take a long time.
    """

    descparts = ["Executable = %s" % script, "Universe = %s" % universe, "Log = %s" % log, "output = %s" % outfile, "error = %s" % errfile, "GetEnv = %s" % True, "Initialdir = %s" % os.getcwd()]

    descparts += ["should_transfer_files = %s" % 'NO', "run_as_owner = %s" % 'True', "Requirements = %s" % requirements, "Notify_user = %s" % mail, "Notification = %s" % 'never', "NiceUser = %s" % niceUser, "Priority = %s" % priority, "Rank = %s" % 'kflops+1000*Memory']

    if arguments:
        descparts.append("Arguments = %s" % arguments)
    # limit the nb of simultaneous runs
    descparts.append("concurrency_limits = %s:%s" % (jobGroup, maxSimultaneousJobsInGroup))
    descparts.append("Queue")

    desc = "\n".join(descparts)
    return submit_text(desc)

def submit_script(script, **kwargs):
    """Like ``submit`` but takes the text of an executable script that
    should be used instead of a filename. Returns the cluster ID along
    with the name of the temporary script file executed. (This should
    probably be removed once the job completes.)
    """
    filename = 'condorpy.jobscript.%s'
    with open(filename, 'w') as f:
        f.write(script)
    os.chmod(filename, 0o755)
    return submit(filename, **kwargs), filename

def wait(jobid, log=LOG_FILE):
    """Waits for a cluster (or specific job) to complete."""
    call("condor_wait %s %s" % (LOG_FILE, str(jobid)))

def getoutput(jobid, log=LOG_FILE, cleanup=True):
    """Waits for a job to complete and then returns its standard output
    and standard error data if the files were given default names.
    Deletes these files after reading them if ``cleanup`` is True.
    """
    wait(jobid, log)
    outfile = OUTFILE_FMT % str(jobid)
    errfile = ERRFILE_FMT % str(jobid)

    with open(outfile) as f:
        stdout = f.read()
    with open(errfile) as f:
        stderr = f.read()

    if cleanup:
        os.unlink(outfile)
        os.unlink(errfile)

    return stdout, stderr

class WaitThread(threading.Thread):
    """A worker that polls Condor log files to observe when jobs
    finish. Each cluster is only waited upon once (after which it is
    "reaped" from the waiting pool).
    """
    def __init__(self, callback, log=LOG_FILE, interval=1):
        """The callable ``callback`` will be invoked with the cluster
        ID of every waited-upon job that finishes. ``interval``
        specifies the polling rate.
        """
        threading.Thread.__init__(self)
        self.callback = callback
        self.log = log
        self.interval = interval
        self.waiting = set()
        self.lock = threading.Lock()
        self.shutdown = False

    def stop(self):
        """Stop the thread soon."""
        with self.lock:
            self.shutdown = True

    def wait(self, clustid):
        """Adds a new job ID to the set of jobs being waited upon."""
        with self.lock:
            self.waiting.add(clustid)

    def run(self):
        while True:
            with self.lock:
                if self.shutdown:
                    return

                # Poll the log file.
                if os.path.exists(self.log):
                    with open(self.log) as f:
                        for line in f:
                            if 'Job terminated.' in line:
                                clustid = re.search(r'\((\d+)\.', line).group(1)
                                clustid = int(clustid)
                                if clustid in self.waiting:
                                    self.callback(clustid)
                                    self.waiting.remove(clustid)

                time.sleep(self.interval)

if __name__ == '__main__':
    ###################
    # Initial example #
    ###################
    # jid = job id (cluster ID of the new job)
    # jfn = job file name
    jid1, jfn1 = submit_script("#!/bin/sh\necho hellow world 1")
    jid2, jfn2 = submit_script("#!/bin/sh\necho hellow world 2")
    try:
        stdout, stderr = getoutput(jid1)
        print "job %i has finished" % jid1
        print stdout
        stdout, stderr = getoutput(jid2)
        print "job %i has finished" % jid2
        print stdout
        print "jobs done"
    finally:
        os.unlink(jfn1)
        os.unlink(jfn2)

    ###################
    ## Second example #
    ###################
    ## myCondor.py should be launched in MagSimus root folder for having good
    ## links
    #listOfJids = []
    #for idxSimu in range(10):
    #    try:
    #        os.mkdir("res/simu1/%s/" % idxSimu)
    #    except:
    #        pass
    #    genesName = 'res/simu1/' + str(idxSimu) + '/genes.%s.list.bz2'
    #    ancGenesName = 'res/simu1/' + str(idxSimu) + '/ancGenes.%s.list.bz2'
    #    jid = submit('src/magSimus1.py',
    #                 arguments='res/speciesTree.phylTree -out:genomeFiles=' + genesName + ' -out:ancGenesFiles=' + ancGenesName + ' -parameterFile=data/parameters.v77 -userRatesFile=data/specRates.v78 +lazyBreakpointAnalyzer',
    #                 universe="vanilla", log=LOG_FILE,
    #                 outfile = OUTFILE_FMT % "$(Cluster)",
    #                 errfile = ERRFILE_FMT % "$(Cluster)")
    #    print "running simu %s (job id %i)" % (idxSimu, jid)
    #    listOfJids.append((idxSimu, jid))

    #for (idxSimu, jid) in listOfJids:
    #    stdout, stderr = getoutput(jid)
    #    # print stderr
    #    print "simu %s done (job id %s)" % (idxSimu, jid)
