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


# It might be a good thing to do execute this two next lines if the ~/condor/ folder is full of old useless files
# cd ~
# rm -rf condor
# mkdir condor

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
from utils import myTools
import stat
import inspect

# Logging messages which are less severe than level will be ignored
# logging.basicConfig(level=logging.DEBUG, format='(%(threadName)-10s) %(message)s')

LOG_FILE = "cdpy.log"
# local buffer folder
USER = getpass.getuser()
REMOTE_BUFF_FOLDER = '/localtmp/' + USER
home = os.path.expanduser("~")
LOCAL_BUFF_FOLDER = home + '/condor'
LASTCOND = LOCAL_BUFF_FOLDER + '/COND'
SCRIPTFILE = "conpy.%s.jobscript"
OUTFILE = "conpy.%s.stdout.log"
ERRFILE = "conpy.%s.stderr.log"
SOE = "conpy.%s.soe"
MACHINES = ["bioclust%02d.bioclust.biologie.ens.fr" % i for i in range(1, 11)] + \
            ["dyoclust%02d.bioclust.biologie.ens.fr" % i for i in range(5, 22)]
           #["dyoclust%02d.bioclust.biologie.ens.fr" % i for i in range(4, 22)]
REQUIREMENT = '(machine != "dyoclust04.bioclust.biologie.ens.fr")'
# Mathieu BAHIN set 10 000 for the total quantity of available simultaneous space for one user.
TOTAL_QUANTITY_AVAILABLE_BY_USER = 10000
MAX_SIMULTANEOUS_JOBS = 100

def quantityEatenByOneJob(maxSimultaneousJobs):
    return float(TOTAL_QUANTITY_AVAILABLE_BY_USER) / float(maxSimultaneousJobs)

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
           requirements='',
           # cf, https://research.cs.wisc.edu/htcondor/CondorWeek2012/presentations/thain-dynamic-slots.pdf
           request_memory='1000',  # request_memory is in mbytes => at least 1Go of RAM
           request_disk='10000',  # request_disk is in kbytes => at least 10Mbytes
           request_cpus='1',
           priority=0,
           # maximum number of simultaneous jobs with the same group name
           maxSimultaneousJobsInGroup=MAX_SIMULTANEOUS_JOBS,
           log=LOG_FILE):
    """Starts a Condor job based on specified parameters. A job
    description is generated. Returns the cluster ID of the new job.

    examples of options:
    1) requirements:
     requirements = (Machine != "bioclust01.bioclust.biologie.ens.fr") && (slotid <= 6)
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
        "request_memory = %s" % request_memory,  # in mbytes
        "request_disk = %s" % request_disk,  # in kbytes
        "request_cpus = %s" % request_cpus,
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


def submit_COND_ManyJobs(condFileForSubmission):
    """Submits a Condor job represented as a job file string. Returns
    the list of jobids of the submitted jobs.
    """
    # print >> sys.stdout, condFileForSubmission
    print >> sys.stderr, "submission of jobs to condor, this may take some time... usually the more jobs the more time (~5sec for 100 jobs and ~25sec for 1000 jobs)"
    out, err = chcall('condor_submit -v', condFileForSubmission)
    print >> sys.stderr, "submission finished"
    # list [..., (clusterId, jid), ...]
    listOfJobids = re.findall(r'Proc (\d+\.\d+)', out)
    for (jobid1, jobid2) in myTools.myIterator.slidingTuple(listOfJobids):
        assert isinstance(jobid1, str) and isinstance(jobid2, str)
        assert jobid1.split('.')[0] == jobid2.split('.')[0]
        assert int(jobid1.split('.')[1]) == int(jobid2.split('.')[1]) - 1
    return listOfJobids

def createCodeFromFunc(pythonFunction, path='./', optimised=True):
        if optimised:
            codeFunc = "#!/usr/bin/python -O\n"
        else:
            codeFunc = "#!/usr/bin/python\n"
        codeFunc += "\n"
        codeFunc += "import sys\n"
        codeFunc += inspect.getsource(pythonFunction)
        codeFunc += "\n"
        codeFunc += "print %s(*sys.argv[1:])\n" % pythonFunction.__name__

        funcNbArgs = inspect.getargspec(pythonFunction)
        funcName = pythonFunction.__name__
        extension = '.py'
        scriptName = funcName + extension
        script = path + scriptName
        with open(script, 'w') as f:
            print >> f, codeFunc
        st = os.stat(script)
        os.chmod(script, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return (script, funcNbArgs)

class CondorSubmitor():
    def __init__(self, executableNameOnCondor):
        assert isinstance(executableNameOnCondor, str)
        self.wrapperExecFileName = LOCAL_BUFF_FOLDER + '/' + executableNameOnCondor + '.sh'
        self.outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE
        self.errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE
        self.signalOfEndFileName = LOCAL_BUFF_FOLDER + '/' + SOE
        self.log = LOG_FILE

    def submit_ManyJobs(self,
                        executable,
                        listOfArguments,
                        universe="vanilla",
                        mail=None,
                        # group name of the job
                        jobGroup='<group>',
                        niceUser=False,
                        requirements=REQUIREMENT,
                        # cf, https://research.cs.wisc.edu/htcondor/CondorWeek2012/presentations/thain-dynamic-slots.pdf
                        request_memory='0.5G',  # request_memory is in mbytes => at least 1Go of RAM
                        request_disk='10000',  # request_disk is in kbytes => at least 10Mbytes
                        request_cpus='1',
                        priority=0,
                        # maximum number of simultaneous jobs with the same group name
                        maxSimultaneousJobsInGroup=MAX_SIMULTANEOUS_JOBS):
        """
        examples of options:
        1) requirements:
        requirements = (Machine != "bioclust01.bioclust.biologie.ens.fr") && (slotid <= 6)
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
            qtEatenByOneJob = 1  # the smallest quantity

        errFileName = self.errFileName % "$(Cluster).%s"
        signalOfEndFileName = self.signalOfEndFileName % "$(Cluster).%s"

        # remove previous log file otherwise the new log will be written at the end and the log file will be long to parse
        for file in [self.log, self.wrapperExecFileName]:
            try:
                os.unlink(file)
            except:
                pass

        # This is useful if the user wrote the executable as a bash command.
        # For instance if the user wrote: 'echo', distutils.spawn.find_executable(executable) returns '/bin/echo'
        # If the user wrote the path toward the executable, this keeps the executable var as a path toward the executable
        executable = distutils.spawn.find_executable(executable)
        with open(self.wrapperExecFileName, 'w') as f:
            print >> f, "#!/bin/bash"
            print >> f, "arrayArgs=(\"$@\")"
            print >> f, "lenArrayArgs=${#arrayArgs[@]}"
            # array contains all the arguments except the two last arguments
            print >> f, "subArray=${arrayArgs[@]:0:$lenArrayArgs-2}"
            print >> f, "%s $subArray > %s" % (executable, self.outFileName % '${arrayArgs[$length-2]}.${arrayArgs[$length-1]}')
            # this line will be sent to the file 'signalOfEndFileName'
            print >> f, "echo \"Signal to inform that job ${arrayArgs[$length-2]}.${arrayArgs[$length-1]} has finished\" "
        # http://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
        st = os.stat(self.wrapperExecFileName)
        os.chmod(self.wrapperExecFileName, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        condFileForSubmission = [
            "Executable = %s" % self.wrapperExecFileName,
            "Universe = %s" % universe,
            "Log = %s" % self.log,
            "GetEnv = %s" % True,
            "Initialdir = %s" % os.getcwd(),
            "should_transfer_files = %s" % 'NO',
            "run_as_owner = %s" % 'True',
            "Requirements = %s" % requirements,
            # cf, https://research.cs.wisc.edu/htcondor/CondorWeek2012/presentations/thain-dynamic-slots.pdf
            "request_cpus = %s" % request_cpus,
            "request_memory = %s" % request_memory,  # in mbytes
            "request_disk = %s" % request_disk,  # in kbytes
            "Notify_user = %s" % mail,
            "Notification = %s" % 'never',
            "NiceUser = %s" % niceUser,
            "Priority = %s" % priority,
            "Rank = %s" % 'kflops+1000*Memory',
            # limit the nb of simultaneous runs
            "concurrency_limits = %s:%s" % (jobGroup, qtEatenByOneJob)
        ]

        for (i, arguments) in enumerate(listOfArguments):
            arguments = arguments + ' ' + '$(Cluster)' + ' ' + str(i)
            condFileForSubmission += "\n"
            if arguments:
                condFileForSubmission += ["Arguments = %s" % arguments]
            condFileForSubmission += ["Input = %s" % '',
                     "Output = %s" % (signalOfEndFileName % i),
                     "Error = %s" % (errFileName % i)]
            condFileForSubmission += ["Queue"]

        condFileForSubmission = "\n".join(condFileForSubmission)
        # write a backup of the last submission file for debug purposes
        with open(LASTCOND, 'w') as f:
            print >> f, condFileForSubmission
        return submit_COND_ManyJobs(condFileForSubmission)

    def waitUntilSignalOfEnd(self, jobid, waitTime=2, maxWaitTime=480):
        signalOfEndFileName = self.signalOfEndFileName % str(jobid)
        while True:
            # if stdout log exists and is non-empty
            if (os.path.isfile(signalOfEndFileName) and os.stat(signalOfEndFileName).st_size != 0) :
                hasFinished = True
                break
            time.sleep(waitTime)
            maxWaitTime -= waitTime
            if maxWaitTime < 0:
                # Case if process timed out
                hasFinished = False
                break
        return (jobid, hasFinished)

    def getoutput_ManyJobs(self, listOfJobids):
        """Waits for a job to complete and then returns its standard output
        and standard error data if the files were given default names.
        """

        pool = Pool()
        for (jobid, hasFinished) in pool.imap_unordered(self.waitUntilSignalOfEnd, tuple(jobid for jobid in listOfJobids)):
            if hasFinished:
                print >> sys.stderr, 'return logs of job', jobid
            else:
                print >> sys.stderr, 'job', jobid, 'has not finished, over max allowed running time'
            signalOfEndFileName = self.signalOfEndFileName % str(jobid)
            try:
                os.remove(signalOfEndFileName)
            except:
                pass
            outFileName = self.outFileName % str(jobid)
            errFileName = self.errFileName % str(jobid)
            yield jobid, outFileName, errFileName

        os.remove(self.wrapperExecFileName)

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

def fib(n):
    n = int(n)
    if n < 2:
        return n
    else:
        return fib(n-2) + fib(n-1)

if __name__ == '__main__':

    ##################
    # first example  #
    ##################

    def fibs_localSequential(nbJobs):
        for n in range(nbJobs):
            print "job%s: fib(35) = %s" % (n, fib(35))

    def fibs_ManyJobs(nbJobs):
        path = './'
        (script, funcNbArgs) = createCodeFromFunc(fib, path=path, optimised=True)
        listOfArguments = [str(35)] * nbJobs
        cs = CondorSubmitor(fib.__name__)
        listOfJobids = cs.submit_ManyJobs(script, listOfArguments, niceUser=True, maxSimultaneousJobsInGroup=None)
        for (jobid, stderrFileName, stdoutFileName) in cs.getoutput_ManyJobs(listOfJobids):
            with open(stdoutFileName, 'r') as f:
                print >> sys.stdout, f.read()
            with open(stderrFileName, 'r') as f:
                print >> sys.stderr, f.read()
            os.unlink(stderrFileName)
            os.unlink(stdoutFileName)
        try:
            os.unlink(script)
        except:
            pass


    ##################
    # second example #
    ##################

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
                       ' -parameterFile=data/parameters.v80 -userRatesFile=data/specRates_MS1.v80 -breakpointAnalyzer'
           listOfArguments.append(arguments)
        cs = CondorSubmitor('magSimus1')
        listOfJids = cs.submit_ManyJobs(executable, listOfArguments, niceUser=True, maxSimultaneousJobsInGroup=None)
        for (jobid, stdoutFileName, stderrFileName) in cs.getoutput_ManyJobs(listOfJids):
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


    import timeit

    nbJobs = 500
    # 1st example
    t_fibs_ManyJobs = timeit.timeit("fibs_ManyJobs(%s)" % nbJobs, setup="from __main__ import fibs_ManyJobs", number=1)
    print >> sys.stderr, "t_fibs_ManyJobs", t_fibs_ManyJobs

    # 2nd example

    # t_magSimus_ManyJobs = timeit.timeit("magSimus_ManyJobs(%s)" % nbJobs, setup="from __main__ import magSimus_ManyJobs", number=1)
    # print >> sys.stderr, "t_magSimus_ManyJobs", t_magSimus_ManyJobs






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

    # over 3000 i get this error:
    #__main__.CommandError: 'condor_submit -v' exited with status 1: '\nWARNING: your Requirements expression refers to TARGET.Memory. This is obsolete. Set request_memory and condor_submit will modify the Requirements expression as needed.\n\nERROR: Failed submission for job 30670.2085 - aborting entire submit\n\nERROR: Failed to queue
    #nbJobs = 3000
    #nbJobs = 200
    # t_magSimus_thread = timeit.timeit("magSimus_thread(%s)" % nbJobs, setup="from __main__ import magSimus_thread", number=1)
    # print >> sys.stderr, "t_magSimus_thread", t_magSimus_thread
    #t_magSimus_buff = timeit.timeit("magSimus_buff(%s)" % nbJobs, setup="from __main__ import magSimus_buff", number=1)
    #print >> sys.stderr, "t_magSimus_buff", t_magSimus_buff

    # print >> sys.stderr, "t_magSimus_ManyJobs", t_magSimus_ManyJobs
    # nbJobs = 200 -> t_magSimus_ManyJobs 60 secs !
    # nbJobs = 1000 -> t_magSimus_ManyJobs 155 secs !
    # nbJobs = 3000 -> t_magSimus_ManyJobs 478 secs !
    # nbJobs = 6000 -> t_magSimus_ManyJobs 914 secs !

    # print >> sys.stderr, "t_magSimus_buff", t_magSimus_buff
    # nbJobs=20 -> t_magSimus_buff 66.1677789688
    # print >> sys.stderr, "t_magSimus_thread", t_magSimus_thread
    # nbJobs=20 -> t_magSimus_thread 50.775026083
    # print >> sys.stderr, "t_magSimus_ManyJobs", t_magSimus_ManyJobs
    # nbJobs=150 -> t_magSimus_thread 45 seconds :D


    #print >> sys.stderr, "Local sequential execution:"
    #t_fibs = timeit.timeit("fibs_localSequential(%s)" % nbJobs, setup="from __main__ import fibs_localSequential", number=1)
    #print >> sys.stderr, "Condor parallel execution:"
    #t_fibs_withThreads = timeit.timeit("fibs_condorThreads(%s)" % nbJobs, setup="from __main__ import fibs_condorThreads", number=1)
    # print >> sys.stderr, "Condor parallel execution (ManyJobs):"


    #print >> sys.stderr, "t_fibs_localSequential", t_fibs
    #print >> sys.stderr, "t_fibs_condorThreads", t_fibs_withThreads

    # nbJobs=500 -> t_fibs_ManyJobs = 40sec
