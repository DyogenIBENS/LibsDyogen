#!/usr/bin/python
# -*- coding: utf-8 -*-

# last author: Joseph Lucas DYOGEN/IBENS 46 rue d'Ulm PARIS 75005
# https://github.com/DyogenIBENS/CondorViaPython

# Fork from
# author:  Adrian Sampson (sampsyo) http://homes.cs.washington.edu/~asampson/
# Found on https://github.com/sampsyo/clusterfutures
# There were no licence so owing to:
#       https://help.github.com/articles/github-terms-of-service/#f
# """... by setting your pages to be viewed publicly, you agree to allow
#    others to view your Content. By setting your repositories to be viewed
#    publicly, you agree to allow others to view and fork your repositories."""
# we forked its deposit https://github.com/DyogenIBENS/CondorViaPython

"""
Abstracts access to a Condor cluster via its command-line tools.
"""
import shutil
import subprocess
import re
import os
import threading
import time
import getpass
import distutils.spawn
import sys
# import logging
from multiprocessing.dummy import Pool


# Logging messages which are less severe than level will be ignored
# logging.basicConfig(level=logging.DEBUG, format='(%(threadName)-10s) %(message)s')

LOG_FILE = "cdpy.log"
# local buffer folder
USER = getpass.getuser()
REMOTE_BUFF_FOLDER = '/localtmp/' + USER
home = os.path.expanduser("~")
LOCAL_BUFF_FOLDER = home + '/condor'
SCRIPTFILE = "conpy.%s.jobscript"
OUTFILE = "conpy.%s.stdout.log"
ERRFILE = "conpy.%s.stderr.log"
MACHINES = ["bioclust%02d.bioclust.biologie.ens.fr" % i for i in range(1, 11)] + \
           ["dyoclust%02d.bioclust.biologie.ens.fr" % i for i in range(4, 22)]
# Mathieu BAHIN set 10 000 for the total quantity of available simultaneous space for one user.
TOTAL_QUANTITY_AVAILABLE_BY_USER = 10000
MAX_SIMULTANEOUS_JOBS = 100

def quantityEatenByOneJob(maxSimultaneousJobs):
    return float(TOTAL_QUANTITY_AVAILABLE_BY_USER) / float(maxSimultaneousJobs)

# TODO: write a method that send jobs via a synthetic COND file with the executable and all info on the top and the
# varying parameters of the different jobs in the following lines

# TODO: use inspection features of python to write a python function into a file
# https://github.com/uqfoundation/dill
# http://stackoverflow.com/questions/1562759/can-python-print-a-function-definition
# def pythonFunctionWrapper(func):
#     local_buff_pythonWrapper = LOCAL_BUFF_FOLDER
#     with open(local_buff_pythonWrapper, 'w') as f:
#         print >> f, "#!/bin/python"
#         print >> f, coreOfFunction
#     os.chmod(local_buff_pythonWrapper, 0o755)

class NotYetFinished(Exception):
    pass

def call(command, stdin=None, returnStderr=True):
    """Invokes a shell command as a subprocess, optionally with some
    data sent to the standard input. Returns the standard output data,
    the standard error, and the return code.

    if returnStderr is None, the sterr value equal None
    """
    if stdin is not None:
        stdin_flag = subprocess.PIPE
    else:
        stdin_flag = None

    if returnStderr:
        stderr_flag = subprocess.PIPE
    else:
        stderr_flag = None

    # FIXME too many open files is due to this line
    proc = subprocess.Popen(command, shell=True, stdin=stdin_flag,
                            stdout=subprocess.PIPE, stderr=stderr_flag, close_fds=True)
    # returns stdout and stderr as strings
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

def wait(jobid, maxTime=None, log=LOG_FILE):
    """Waits for a cluster (or specific job) to complete.
    maxTime: Wait no more than this time (None means unlimited)"""
    if maxTime:
        (stdout, _, _) = call("condor_wait -wait %s %s %s" % (maxTime, log, str(jobid)))
        # logging.debug('in wait' + stdout + '. Is it "Time expired." or "All jobs done."')
        if stdout == 'Time expired.':
            finished = False
        else:
            assert stdout == 'All jobs done.'
            finished = True
    else:
        call("condor_wait %s %s" % (log, str(jobid)))
        finished = True
    return finished

def wait_Jid((jobid, maxTime, log)):
    """Waits for a cluster (or specific job) to complete.
    maxTime: Wait no more than this time (None means unlimited)"""
    if maxTime:
        (stdout, _, _) = call("condor_wait -wait %s %s %s" % (maxTime, log, str(jobid)))
        # logging.debug('in wait' + stdout + '. Is it "Time expired." or "All jobs done."')
        if stdout == 'Time expired.':
            finished = False
        else:
            assert stdout == 'All jobs done.'
            finished = True
    else:
        call("condor_wait %s %s" % (log, str(jobid)))
        finished = True
    return jobid

def hasFinished(jobid, log=LOG_FILE):
    """Return the immediately the current job status"""
    (stdout, _, _) = call("condor_wait -wait %s %s %s" % (0, log, str(jobid)))
    if stdout[-1] == 'All jobs done':
        finished = True
    else:
        finished = False
    return finished

def submit_COND_OneJob(COND):
    """Submits a Condor job represented as a job file string. Returns
    the cluster ID of the submitted job.
    """
    out, _ = chcall('condor_submit -v', COND)
    jobid = re.search(r'Proc (\d+)\.0', out).group(1)
    return int(jobid)

def submit_OneJob(executable, universe="vanilla",
           arguments=None,
           mail=None,
           # group name of the job
           jobGroup='<group>',
           niceUser=False,
           requirements='(Memory > 1024)',
           priority=0,
           # maximum number of simultaneous jobs with the same group name
           maxSimultaneousJobsInGroup=MAX_SIMULTANEOUS_JOBS,
           log=LOG_FILE):
    """Starts a Condor job based on specified parameters. A job
    description is generated. Returns the cluster ID of the new job.

    examples of options:
    1) requirements:
    requirements = (Memory > 1024) && (Machine != "bioclust01.bioclust.biologie.ens.fr") && (slotid <= 6)
    (slotid <= 6) means that the job won't occupy more than 6 thread on a multi-core machine.
    This option may be interesting if the job uses several threads on a same machine.
    For instance Blast uses 4 threads and it would be interesting to send the job with (slotid <= 4) on an octo-core
    machine.
    2) priority:
    Each job may have a priority level specified by 'priority', an int >= 0.
    Condor executes jobs by decreasing priority levels.
    If a job1 has a priority1 and job2 has a priority2 and if priority1 > priority2, job1 will be executed before job2
    This won't change the priority between two different users, it just changes the order of jobs of a same user.
    3) niceUser:
    If the user is planning to send a big amount of jobs that do not take long to execute, he may execute them with
    the niceUser option set to True (and no limit of simultaneous runs), thus all other users will have the priority
    over him.
    4) maxSimRuns is an upper limit of simultaneous runs on all the machines of the cluster. The user can send an
    unlimited amount of jobs that will be queued but no more than 100 jobs will be executed simultaneously.
    This limit limit the saturation of cores and keep some cores free for other users.
    This option is especially important if jobs take a long time.
    """

    if maxSimultaneousJobsInGroup is not None:
        qtEatenByOneJob = quantityEatenByOneJob(maxSimultaneousJobsInGroup)
    else:
        qtEatenByOneJob = 1 # the smallest quantity

    # This is useful if the user wrote the executable as a bash command.
    # For instance if the user wrote: 'echo', distutils.spawn.find_executable(executable) returns '/bin/echo'
    # If the user wrote the path toward the executable, this keeps the executable var as a path toward the executable
    executable = distutils.spawn.find_executable(executable)

    outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % "$(Cluster)"
    errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % "$(Cluster)"

    COND = [
        "Executable = %s" % executable,
        "Universe = %s" % universe,
        "Log = %s" % log,
        "GetEnv = %s" % True,
        "Initialdir = %s" % os.getcwd(),
        "should_transfer_files = %s" % 'NO',
        "run_as_owner = %s" % 'True',
        "Requirements = %s" % requirements,
        "Notify_user = %s" % mail,
        "Notification = %s" % 'never',
        "NiceUser = %s" % niceUser,
        "Priority = %s" % priority,
        "Rank = %s" % 'kflops+1000*Memory',
        # limit the nb of simultaneous runs
        "concurrency_limits = %s:%s" % (jobGroup, qtEatenByOneJob)
    ]

    COND += ["\n"]

    if arguments:
        COND += ["Arguments = %s" % arguments]
    COND += ["Input = %s" % '',
             "Output = %s" % outFileName,
             "Error = %s" % errFileName]

    COND += ["Queue"]

    COND = "\n".join(COND)
    return submit_COND_OneJob(COND)


def getoutput(jobid, waitMaxTime=None, log=LOG_FILE):
    """Waits for a job to complete and then returns its standard output
    and standard error data if the files were given default names.
    Deletes these files after reading them if ``cleanup`` is True.
    """

    outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % str(jobid)
    errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % str(jobid)

    finished = wait(jobid, maxTime=waitMaxTime, log=log)

    if not finished:
        raise NotYetFinished

    return outFileName, errFileName


def submitWithBuffer(command, jobName, **kwargs):
    try:
        # create a local buff dir
        os.mkdir(LOCAL_BUFF_FOLDER)
    except OSError as e:
        assert e.strerror == 'File exists'
    except Exception as e:
        raise e

    local_buff_script = LOCAL_BUFF_FOLDER + '/' + SCRIPTFILE % jobName
    remote_buff_outfile = REMOTE_BUFF_FOLDER + '/' + OUTFILE % jobName
    remote_buff_errfile = REMOTE_BUFF_FOLDER + '/' + ERRFILE % jobName
    local_buff_outfile = LOCAL_BUFF_FOLDER + '/' + OUTFILE % jobName
    local_buff_errfile = LOCAL_BUFF_FOLDER + '/' + ERRFILE % jobName

    with open(local_buff_script, 'w') as f:
        # prepare the script that will make the dir on the remote cluster machine
        print >> f, "#!/bin/sh"
        print >> f, "mkdir %s 2>/dev/null" % REMOTE_BUFF_FOLDER
        print >> f, command + " > " + remote_buff_outfile + " 2> " + remote_buff_errfile
        # Once the script is finished, move the results in the wished directory
        print >> f, "mv " + remote_buff_outfile + " " + local_buff_outfile
        print >> f, "mv " + remote_buff_errfile + " " + local_buff_errfile
    os.chmod(local_buff_script, 0o755)

    jid = submit_OneJob(local_buff_script, **kwargs)

    return jid


def getOutputWithBuffer(jobid, jobName, waitMaxTime=None, log=LOG_FILE, cleanup=True):

    local_buff_script = LOCAL_BUFF_FOLDER + '/' + SCRIPTFILE % jobName

    local_buff_outfile_jid = LOCAL_BUFF_FOLDER + '/' + OUTFILE % str(jobid)
    local_buff_errfile_jid = LOCAL_BUFF_FOLDER + '/' + ERRFILE % str(jobid)
    outfileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % jobName
    errfileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % jobName

    finished = wait(jobid, maxTime=waitMaxTime, log=log)

    if not finished:
        raise NotYetFinished

    if cleanup:
        # remove temporary files
        os.unlink(local_buff_script)

        os.unlink(local_buff_outfile_jid)
        os.unlink(local_buff_errfile_jid)

    return (outfileName, errfileName)


def submit_COND_ManyJobs(COND):
    """Submits a Condor job represented as a job file string. Returns
    the list of jobids of the submitted jobs.
    """
    # print >> sys.stdout, COND
    print >> sys.stderr, "submission of jobs to condor, this may take some time... usually the more jobs the more time (~5sec for 100 jobs and ~25sec for 1000 jobs)"
    out, err = chcall('condor_submit -v', COND)
    print >> sys.stderr, "submission finished"
    # list [..., (clusterId, jid), ...]
    listOfJobids = re.findall(r'Proc (\d+\.\d+)', out)
    return listOfJobids


def submit_ManyJobs(executable,
                    listOfArguments,
                    universe="vanilla",
                    mail=None,
                    # group name of the job
                    jobGroup='<group>',
                    niceUser=False,
                    requirements='(Memory > 1024)',
                    priority=0,
                    # maximum number of simultaneous jobs with the same group name
                    maxSimultaneousJobsInGroup=MAX_SIMULTANEOUS_JOBS,
                    log=LOG_FILE):
    """Starts a Condor job based on specified parameters. A job
    description is generated. Returns the cluster ID of the new job.

    examples of options:
    1) requirements:
    requirements = (Memory > 1024) && (Machine != "bioclust01.bioclust.biologie.ens.fr") && (slotid <= 6)
    (slotid <= 6) means that the job won't occupy more than 6 thread on a multi-core machine.
    This option may be interesting if the job uses several threads on a same machine.
    For instance Blast uses 4 threads and it would be interesting to send the job with (slotid <= 4) on an octo-core
    machine.
    2) priority:
    Each job may have a priority level specified by 'priority', an int >= 0.
    Condor executes jobs by decreasing priority levels.
    If a job1 has a priority1 and job2 has a priority2 and if priority1 > priority2, job1 will be executed before job2
    This won't change the priority between two different users, it just changes the order of jobs of a same user.
    3) niceUser:
    If the user is planning to send a big amount of jobs that do not take long to execute, he may execute them with
    the niceUser option set to True (and no limit of simultaneous runs), thus all other users will have the priority
    over him.
    4) maxSimRuns is an upper limit of simultaneous runs on all the machines of the cluster. The user can send an
    unlimited amount of jobs that will be queued but no more than 100 jobs will be executed simultaneously.
    This limit limit the saturation of cores and keep some cores free for other users.
    This option is especially important if jobs take a long time.
    """

    if maxSimultaneousJobsInGroup is not None:
        qtEatenByOneJob = quantityEatenByOneJob(maxSimultaneousJobsInGroup)
    else:
        qtEatenByOneJob = 1 # the smallest quantity

    # This is useful if the user wrote the executable as a bash command.
    # For instance if the user wrote: 'echo', distutils.spawn.find_executable(executable) returns '/bin/echo'
    # If the user wrote the path toward the executable, this keeps the executable var as a path toward the executable

    # remove previous log file otherwise the new log will be written at the end and the log file will be long to parse
    outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % "$(Cluster).%s"
    errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % "$(Cluster).%s"
    try:
        os.unlink(log)
    except:
        pass
    try:
        os.unlink(outFileName)
    except:
        pass
    try:
        os.unlink(errFileName)
    except:
        pass

    executable = distutils.spawn.find_executable(executable)

    COND = [
        "Executable = %s" % executable,
        "Universe = %s" % universe,
        "Log = %s" % log,
        "GetEnv = %s" % True,
        "Initialdir = %s" % os.getcwd(),
        "should_transfer_files = %s" % 'NO',
        "run_as_owner = %s" % 'True',
        "Requirements = %s" % requirements,
        "Notify_user = %s" % mail,
        "Notification = %s" % 'never',
        "NiceUser = %s" % niceUser,
        "Priority = %s" % priority,
        "Rank = %s" % 'kflops+1000*Memory',
        # limit the nb of simultaneous runs
        "concurrency_limits = %s:%s" % (jobGroup, qtEatenByOneJob)
    ]

    for (i, arguments) in enumerate(listOfArguments):
        COND += "\n"
        if arguments:
            COND += ["Arguments = %s" % arguments]
        COND += ["Input = %s" % '',
                 "Output = %s" % (outFileName % i),
                 "Error = %s" % (errFileName % i)]
        COND += ["Queue"]

    COND = "\n".join(COND)
    return submit_COND_ManyJobs(COND)

def waitUntilLogOutExists(jobid, waitTime=2):
    outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % str(jobid)
    while True:
        # if the file exists and is not empty
        if os.path.isfile(outFileName) and os.stat(outFileName).st_size != 0:
            return jobid
        time.sleep(waitTime)

def getoutput_ManyJobs(listOfJobids):
    """Waits for a job to complete and then returns its standard output
    and standard error data if the files were given default names.
    Deletes these files after reading them if ``cleanup`` is True.
    """

    pool = Pool()
    for jobid in pool.imap_unordered(waitUntilLogOutExists, (jobid for jobid in listOfJobids)):
        print >> sys.stderr, 'return logs of job', jobid
        outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % str(jobid)
        errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % str(jobid)
        yield jobid, outFileName, errFileName

def execBashCmdOnAllMachines(command, machines=MACHINES, mail=None, log=LOG_FILE):
    """
    If you have pbs on your jobs, the directory REMOTE_BUFF_FOLDER won't be clean automatically.
    this script allows to remove all files in REMOTE_BUFF_FOLDER on all machines
    """
    listOfJids = []
    for (i, machine) in enumerate(machines):
        jobName = str(i)
        jid = submitWithBuffer(command, jobName,
                               mail=mail,
                               requirements="(Machine==\"%s\")" % machine,
                               log=log)
        listOfJids.append((machine, jid, jobName))

    for (machine, jid, jobName) in listOfJids:
        print >> sys.stderr, (machine, jid, jobName)
        stdoutFileName, stderrFileName = getOutputWithBuffer(jid, jobName)
        printFileIntoStream(stdoutFileName, sys.stdout)
        printFileIntoStream(stderrFileName, sys.stderr)
        os.unlink(stdoutFileName)
        os.unlink(stderrFileName)

# log_stderr = logging.getLogger('stderr')
# log_stderr.propagate = False
# log_stderr.setLevel(logging.DEBUG)
# #log_stderr.addHandler(logging.StreamHandler(stream=sys.stderr))
#
# log_stdout = logging.getLogger('stdout')
# log_stdout.propagate = False
# log_stdout.setLevel(level=logging.INFO)
#log_stdout.addHandler(logging.StreamHandler(stream=sys.stdout))

# user may want to log flux to files
# # add a file handler
# fh = logging.FileHandler('myapp.log')
# fh.setLevel(logging.WARNING)
# # create a formatter and set the formatter for the handler.
# frmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# fh.setFormatter(frmt)
# # add the Handler to the logger
# lgr.addHandler(fh)

#  """
# https://docs.python.org/2/library/logging.html
# The logging module is intended to be thread-safe without any special work needing to be done by its clients.
# It achieves this though using threading locks; there is one lock to serialize access to the module’s shared data,
# and each handler also creates a lock to serialize access to its underlying I/O.
#
# If you are implementing asynchronous signal handlers using the signal module, you may not be able to use logging
# from within such handlers. This is because lock implementations in the threading module are not always re-entrant,
# and so cannot be invoked from such signal handlers.
#      """

class condorThread(threading.Thread):
    def __init__(self, command=None, monitoringTimeStep=5, name=None, callBack=None, log=LOG_FILE, verbose=False):
        """*target* is the callable object to be invoked by the run()
        method. Defaults to None, meaning nothing is called.

        *name* is the thread name. By default, a unique name is constructed of
        the form "Thread-N" where N is a small decimal number.

        *args* is the argument tuple for the target invocation. Defaults to ().

        *kwargs* is a dictionary of keyword arguments for the target
        invocation. Defaults to {}."""
        assert isinstance(command, str)
        self.command = command
        self.target = command.split(' ')[0]
        self.args = tuple(command.split(' ')[1:])
        self.monitoringTimeStep = monitoringTimeStep
        """A string used for identification purposes only. It has no semantics. Multiple threads may be given the same
           name. The initial name is set by the constructor.
        """
        self.tName = name
        self.log = log
        self.verbose = verbose
        self.res = None
        self.callBack = callBack
        threading.Thread.__init__(self, group=None, target=self.target, name=self.tName,
                                  args=self.args, kwargs=None)

    # def logErr(self, msg):
    #     if self.verbose:
    #         log_stderr.log(logging.DEBUG, msg)
    #     else:
    #         pass
    #
    # def logOut(self, msg):
    #     log_stdout.log(logging.INFO, msg)

    # def logErr(self, msg):
    #     if self.verbose:
    #         logging.debug(msg)

    def run(self):
        # self.logErr('run')
        if self.command:
            self.jid = submitWithBuffer(self.command, self.tName,
                                        mail=None,
                                        # group name of the job
                                        jobGroup='<group>',
                                        niceUser=False,
                                        requirements='(Memory > 1024)',
                                        priority=0,
                                        # maximum number of simultaneous jobs with the same group name
                                        maxSimultaneousJobsInGroup=quantityEatenByOneJob(MAX_SIMULTANEOUS_JOBS),
                                        log=self.log)
            while True:
                # self.logErr('Monitoring')
                time.sleep(self.monitoringTimeStep)
                try:
                    (self.outFileName, self.errFileName) = getOutputWithBuffer(self.jid, self.tName, waitMaxTime=0, log=self.log)
                    self.res = (self.outFileName, self.errFileName)
                    # self.logErr("thread %s has finished" % self.tName)
                    if self.callBack:
                        self.callBack(self.res)
                    break
                except NotYetFinished:
                    # self.logErr('NotYetFinished')
                    continue
                except Exception as e:
                    # self.logErr(str(type(e)) + ' ' + str(e))
                    raise e

    def __del__(self):
        del self.command, self.target, self.args, self.log, self.res

def printLogErr((stdoutFileName, stderrFileName)):
    printFileIntoStream(stderrFileName, sys.stderr)

def printLogOut((stdoutFileName, stderrFileName)):
    printFileIntoStream(stdoutFileName, sys.stdout)

def printLogs((stdoutFileName, stderrFileName)):
    printLogErr((stdoutFileName, stdoutFileName))
    printLogOut((stderrFileName, stderrFileName))

def printFileIntoStream(fileName, stream):
    with open(fileName, 'r') as f:
        try:
            print >> stream, f.read()[:-1]
        except:
            pass

if __name__ == '__main__':

    ##################
    #  first example #
    ##################
    def helloWorld_noBuff(nbJobs):
        listOfJNamesIds = []
        for jobName in range(nbJobs):
            # jid = job id (cluster ID of the new job)
            jid = submit_OneJob("echo", arguments="Hello job %s!" % jobName)
            print "job %s has id %i" % (jobName, jid)
            listOfJNamesIds.append((jobName, jid))

        for (jobName, jid) in listOfJNamesIds:
            (stdoutFileName, stderrFileName) = getoutput(jid)
            print "job %s has finished" % jobName
            printFileIntoStream(stdoutFileName, sys.stdout)
            # printFileIntoStream(stderrFileName, sys.stderr)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)
        print "jobs done"

    def helloWorld_buff(nbJobs):
        # Same with buffer
        listOfJNamesIds = []
        for jobName in range(nbJobs):
            jid = submitWithBuffer("echo \"Hello job %s!\"" % jobName, jobName)
            print "job %s has id %i" % (jobName, jid)
            listOfJNamesIds.append((jobName, jid))

        for (jobName, jid) in listOfJNamesIds:
            (stdoutFileName, stderrFileName) = getOutputWithBuffer(jid, jobName)
            print "job %s has finished" % jobName
            printFileIntoStream(stdoutFileName, sys.stdout)
            # printFileIntoStream(stderrFileName, sys.stderr)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)
        print "jobs done"

    def helloWorld_thread(nbJobs):
        # Same with threads
        verbose = True
        listOfThreads = []
        for jobName in range(nbJobs):
            t = condorThread(command="echo \"Hello job %s!\"" % jobName, name=jobName, callBack=printLogOut, verbose=verbose)
            t.start()
            listOfThreads.append(t)

        for t in listOfThreads:
            t.join()
            (stdoutFileName, stderrFileName) = t.res
            # printFileIntoStream(stdoutFileName, sys.stdout)
            # printFileIntoStream(stderrFileName, sys.stderr)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)
            del t

    # ##################
    # # second example #
    # ##################
    def magSimus_buff(nbJobs):
        # This example uses a software called src/magSimus1.py but with specific
        # arguments but it could use any other executable and arguments.
        # condor.py should be launched in MagSimus root folder for having good links
        listOfJids = []
        for idxSimu in range(nbJobs):
           try:
               os.mkdir("res/simu1/%s/" % idxSimu)
           except:
               pass
           jobName = idxSimu
           genesName = 'res/simu1/' + str(idxSimu) + '/genes.%s.list.bz2'
           ancGenesName = 'res/simu1/' + str(idxSimu) + '/ancGenes.%s.list.bz2'
           executable = 'src/magSimus1.py'
           arguments = 'res/speciesTree.phylTree -out:genomeFiles=' + genesName +\
                       ' -out:ancGenesFiles=' + ancGenesName +\
                       ' -parameterFile=data/parameters.v80 -userRatesFile=data/specRates_MS1.v80 +lazyBreakpointAnalyzer'
           command = executable + ' ' + arguments
           jid = submitWithBuffer(command, jobName)
           print "simu %s (job id %i) sent to condor" % (idxSimu, jid)
           listOfJids.append((jid, jobName))

        for (jid, jobName) in listOfJids:
            idxSimu = jobName
            stdoutFileName, stderrFileName = getOutputWithBuffer(jid, jobName)
            print "data of %s returned (thread %s)" % (idxSimu, jobName)
            # printFileIntoStream(stdoutFileName, sys.stdout)
            printFileIntoStream(stderrFileName, sys.stderr)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)

    def magSimus_thread(nbJobs):
        # Same with threads
        listOfThreads = []
        for idxSimu in range(nbJobs):
           try:
               os.mkdir("res/simu1/%s/" % idxSimu)
           except:
               pass
           jobName = 'simu' + str(idxSimu)
           genesName = 'res/simu1/' + str(idxSimu) + '/genes.%s.list.bz2'
           ancGenesName = 'res/simu1/' + str(idxSimu) + '/ancGenes.%s.list.bz2'
           executable = 'src/magSimus1.py'
           arguments = 'res/speciesTree.phylTree -out:genomeFiles=' + genesName +\
                       ' -out:ancGenesFiles=' + ancGenesName +\
                       ' -parameterFile=data/parameters.v80 -userRatesFile=data/specRates_MS1.v80 +lazyBreakpointAnalyzer'
           command = executable + ' ' + arguments

           # do not print stdout because here we do not need it
           t = condorThread(command=command, name=jobName, callBack=printLogErr, verbose=True)
           t.start()
           print "simu %s (thread name %s) sent to condor" % (idxSimu, t.tName)
           listOfThreads.append((t, idxSimu))

        for (t, idxSimu) in listOfThreads:
            t.join()
            stdoutFileName, stderrFileName = t.res
            print "data of %s returned (thread %s)" % (idxSimu, t.tName)
            # printFileIntoStream(stdoutFileName, sys.stdout)
            # The callBack function already print the errLog
            # printFileIntoStream(stderrFileName, sys.stderr)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)
            del t

    def magSimus_ManyJobs(nbJobs):
        # Same with threads
        executable = 'src/magSimus1.py'
        listOfArguments = []
        for idxSimu in range(nbJobs):
           try:
               os.mkdir("res/simu1/%s/" % idxSimu)
           except:
               pass
           genesName = 'res/simu1/' + str(idxSimu) + '/genes.%s.list.bz2'
           ancGenesName = 'res/simu1/' + str(idxSimu) + '/ancGenes.%s.list.bz2'
           arguments = 'res/speciesTree.phylTree -out:genomeFiles=' + genesName +\
                       ' -out:ancGenesFiles=' + ancGenesName +\
                       ' -parameterFile=data/parameters.v80 -userRatesFile=data/specRates_MS1.v80 +lazyBreakpointAnalyzer'
           listOfArguments.append(arguments)
        listOfJids = submit_ManyJobs(executable, listOfArguments, niceUser=True, maxSimultaneousJobsInGroup=None)
        for (jobid, stdoutFileName, stderrFileName) in getoutput_ManyJobs(listOfJids):
            # print the 3 first lines of stdout and stderr logs
            with open(stdoutFileName, 'r') as f:
                print >> sys.stdout, f.readline(),
                print >> sys.stdout, f.readline(),
                print >> sys.stdout, f.readline(),
            with open(stderrFileName, 'r') as f:
                print >> sys.stderr, f.readline(),
                print >> sys.stdout, f.readline(),
                print >> sys.stdout, f.readline(),
            idxSimu = jobid.split('.')[1]
            shutil.move(stdoutFileName, 'res/simu1/' + str(idxSimu) +'/logOut')
            shutil.move(stderrFileName, 'res/simu1/' + str(idxSimu) +'/logErr')
            # os.unlink(stdoutFileName)
            # os.unlink(stderrFileName)

    #######################################
    # Execute commands on remote machines #
    #######################################
    def example3():
        # echo 'toto
        execBashCmdOnAllMachines('echo "toto"')

        # create a file toto in remote folders
        execBashCmdOnAllMachines('touch ' + REMOTE_BUFF_FOLDER + '/toto')
        # look inside the remote temporary folders
        # To check if all temporary files have been removed properly
        execBashCmdOnAllMachines('ls ' + REMOTE_BUFF_FOLDER)

        # remove all files in remote folders
        execBashCmdOnAllMachines('rm ' + REMOTE_BUFF_FOLDER + '/*')

    def createFib35Script():
            code =\
            ("#!/usr/bin/python\n"
             "def fib(n):\n"
             "    if n < 2:\n"
             "        return n\n"
             "    return fib(n-2) + fib(n-1)\n"
             "print 'fib(35) =', fib(35)"
             )
            filename = './fib35.py'
            try:
                os.unlink(filename)
            except:
                pass
            with open(filename, 'w') as f:
                print >> f, code
            os.chmod(filename, 0755)

    def fibs_localSequential(nbJobs):
        def fib(n):
            if n < 2:
                return n
            return fib(n-2) + fib(n-1)

        for n in range(nbJobs):
            print "job%s: fib(35) = %s" % (n, fib(35))

    def fibs_condorThreads(nbJobs):

        listOfThreads = []
        createFib35Script()
        for n in range(nbJobs):
           jobName = 'fib35_' + str(n)
           command = 'python ./fib35.py'
           t = condorThread(command=command, name=jobName, verbose=False)
           t.start()
           assert t.tName == jobName
           listOfThreads.append((t, jobName))

        for (t, jobName) in listOfThreads:
            t.join()
            stdoutFileName, stderrFileName = t.res
            print >> sys.stdout, jobName + ':',
            printFileIntoStream(stdoutFileName, sys.stdout)
            # remove files
            os.unlink(stdoutFileName)
            os.unlink(stderrFileName)
            try:
                os.unlink('./fib35.py')
            except:
                pass
            del t

    def fibs_ManyJobs(nbJobs):
        createFib35Script()

        command = './fib35.py'
        listOfArguments = [None] * nbJobs
        listOfJobids = submit_ManyJobs(command, listOfArguments, niceUser=True, maxSimultaneousJobsInGroup=None)
        for (jobid, stderrFileName, stdoutFileName) in getoutput_ManyJobs(listOfJobids):
            with open(stdoutFileName, 'r') as f:
                print >> sys.stdout, f.read()
            with open(stderrFileName, 'r') as f:
                print >> sys.stderr, f.read()
            os.unlink(stderrFileName)
            os.unlink(stdoutFileName)
        try:
            os.unlink('./fib35.py')
        except:
            pass


    import timeit

    #nbJobs = 20
    #t_helloWorld_noBuff = timeit.timeit("helloWorld_noBuff(%s)" % nbJobs, setup="from __main__ import helloWorld_noBuff", number=1)
    #print >> sys.stderr, "helloWorld_noBuff", t_helloWorld_noBuff
    #t_helloWorld_buff = timeit.timeit("helloWorld_buff(%s)" % nbJobs, setup="from __main__ import helloWorld_buff", number=1)
    #print >> sys.stderr, "helloWorld_buff", t_helloWorld_buff
    #t_helloWorld_thread = timeit.timeit("helloWorld_thread(%s)" % nbJobs, setup="from __main__ import helloWorld_thread", number=1)
    #print >> sys.stderr, "t_helloWorld_thread", t_helloWorld_thread

    #print >> sys.stderr, "helloWorld_noBuff", t_helloWorld_noBuff
    ## helloWorld_noBuff 28.5351910591
    #print >> sys.stderr, "helloWorld_buff", t_helloWorld_buff
    ## helloWorld_buff 46.931251049
    #print >> sys.stderr, "t_helloWorld_thread", t_helloWorld_thread
    ## t_helloWorld_thread 12.5498409271

    nbJobs = 3000
    #nbJobs = 200
    # t_magSimus_thread = timeit.timeit("magSimus_thread(%s)" % nbJobs, setup="from __main__ import magSimus_thread", number=1)
    # print >> sys.stderr, "t_magSimus_thread", t_magSimus_thread
    #t_magSimus_buff = timeit.timeit("magSimus_buff(%s)" % nbJobs, setup="from __main__ import magSimus_buff", number=1)
    #print >> sys.stderr, "t_magSimus_buff", t_magSimus_buff
    t_magSimus_ManyJobs = timeit.timeit("magSimus_ManyJobs(%s)" % nbJobs, setup="from __main__ import magSimus_ManyJobs", number=1)
    print >> sys.stderr, "t_magSimus_ManyJobs", t_magSimus_ManyJobs
    # nbJobs = 1000 -> t_magSimus_ManyJobs 155 secs !
    # nbJobs = 3000 -> t_magSimus_ManyJobs 478 secs !

    # print >> sys.stderr, "t_magSimus_buff", t_magSimus_buff
    # nbJobs=20 -> t_magSimus_buff 66.1677789688
    # print >> sys.stderr, "t_magSimus_thread", t_magSimus_thread
    # nbJobs=20 -> t_magSimus_thread 50.775026083
    # print >> sys.stderr, "t_magSimus_ManyJobs", t_magSimus_ManyJobs
    # nbJobs=150 -> t_magSimus_thread 45 seconds :D

    # nbJobs = 500
    #print >> sys.stderr, "Local sequential execution:"
    #t_fibs = timeit.timeit("fibs_localSequential(%s)" % nbJobs, setup="from __main__ import fibs_localSequential", number=1)
    #print >> sys.stderr, "Condor parallel execution:"
    #t_fibs_withThreads = timeit.timeit("fibs_condorThreads(%s)" % nbJobs, setup="from __main__ import fibs_condorThreads", number=1)
    # print >> sys.stderr, "Condor parallel execution (ManyJobs):"
    # t_fibs_ManyJobs = timeit.timeit("fibs_ManyJobs(%s)" % nbJobs, setup="from __main__ import fibs_ManyJobs", number=1)

    #print >> sys.stderr, "t_fibs_localSequential", t_fibs
    #print >> sys.stderr, "t_fibs_condorThreads", t_fibs_withThreads
    # print >> sys.stderr, "t_fibs_ManyJobs", t_fibs_ManyJobs
    # nbJobs=500 -> t_fibs_ManyJobs = 40sec
