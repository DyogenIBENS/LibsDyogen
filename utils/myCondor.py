#!/usr/bin/python
# -*- coding: utf-8 -*-

# last author: Joseph Lucas DYOGEN/IBENS 46 rue d'Ulm PARIS 755
# https://github.com/DyogenIBENS/CondorViaPython

# Fork from
# author:  Adrian Sampson (sampsyo) http://homes.cs.washington.edu/~asampson/
# Found on https://github.com/sampsyo/clusterfutures
# There were no licence so owing to:
#       https://help.github.com/articles/github-terms-of-service/#f
# """... by setting your pages to be viewed publicly, you agree to allow
#    others to view your Content. By setting your repositories to be viewed
#    publicly, you agree to allow others to view and fork your repositories."""
# we forked its deposit https://github.com/andeElliott/condorPythonShell

"""
Abstracts access to a Condor cluster via its command-line tools.
"""

import subprocess
import re
import os
import threading
import time
import getpass
import distutils.spawn
import sys
import logging


# Logging messages which are less severe than level will be ignored
logging.basicConfig(level=logging.DEBUG, format='(%(threadName)-10s) %(message)s')

LOG_FILE = "condorpy.log"
# local buffer folder
USER = getpass.getuser()
REMOTE_BUFF_FOLDER = '/localtmp/' + USER
home = os.path.expanduser("~")
LOCAL_BUFF_FOLDER = home + '/condor'
SCRIPTFILE = "condorpy.%s.jobscript"
OUTFILE = "condorpy.%s.stdout.log"
ERRFILE = "condorpy.%s.stderr.log"
# dyoclust07-08 are not 100% functional
MACHINES = ["bioclust%02d.bioclust.biologie.ens.fr" % i for i in range(1, 11)] + \
           ["dyoclust%02d.bioclust.biologie.ens.fr" % i for i in range(4, 7) + range(9, 22)]


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

def wait(jobid, maxTime=None, log=LOG_FILE):
    """Waits for a cluster (or specific job) to complete.
    maxTime: Wait no more than this time (None means unlimited)"""
    if maxTime:
        (stdout, _, _) = call("condor_wait -wait %s %s %s" % (maxTime, log, str(jobid)))
        logging.debug('in wait' + stdout + '. Is it "Time expired." or "All jobs done."')
        if stdout == 'Time expired.':
            finished = False
        else:
            assert stdout == 'All jobs done.'
            finished = True
    else:
        call("condor_wait %s %s" % (log, str(jobid)))
        finished = True
    return finished

def hasFinished(jobid, log=LOG_FILE):
    """Return the immediately the current job status"""
    (stdout, _, _) = call("condor_wait -wait %s %s %s" % (0, log, str(jobid)))
    if stdout[-1] == 'All jobs done':
        finished = True
    else:
        finished = False
    return finished

def submit_text(job):
    """Submits a Condor job represented as a job file string. Returns
    the cluster ID of the submitted job.
    """
    out, _ = chcall('condor_submit -v', job)
    jobid = re.search(r'Proc (\d+)\.0', out).group(1)
    return int(jobid)

def submit(executable, universe="vanilla",
           arguments=None,
           mail=None,
           # group name of the job
           jobGroup='<group>',
           niceUser=False,
           requirements='(Memory > 1024)',
           priority=0,
           # maximum number of simultaneous jobs with the same group name
           maxSimultaneousJobsInGroup=100,
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

    # This is useful if the user wrote the executable as a bash command.
    # For instance if the user wrote: 'echo', distutils.spawn.find_executable(executable) returns '/bin/echo'
    # If the user wrote the path toward the executable, this keeps the executable var as a path toward the executable
    executable = distutils.spawn.find_executable(executable)

    outFileName = LOCAL_BUFF_FOLDER + '/' + OUTFILE % "$(Cluster)"
    errFileName = LOCAL_BUFF_FOLDER + '/' + ERRFILE % "$(Cluster)"

    descparts = [
        "Executable = %s" % executable,
        "Universe = %s" % universe,
        "Log = %s" % log,
        "output = %s" % outFileName,
        "error = %s" % errFileName,
        "GetEnv = %s" % True,
        "Initialdir = %s" % os.getcwd()
    ]

    descparts += [
        "should_transfer_files = %s" % 'NO',
        "run_as_owner = %s" % 'True',
        "Requirements = %s" % requirements,
        "Notify_user = %s" % mail,
        "Notification = %s" % 'never',
        "NiceUser = %s" % niceUser,
        "Priority = %s" % priority,
        "Rank = %s" % 'kflops+1000*Memory'
    ]

    if arguments:
        descparts.append("Arguments = %s" % arguments)
    # limit the nb of simultaneous runs
    descparts.append("concurrency_limits = %s:%s" % (jobGroup, maxSimultaneousJobsInGroup))
    descparts.append("Queue")

    desc = "\n".join(descparts)
    return submit_text(desc)

def getoutput(jobid, waitMaxTime=None, log=LOG_FILE, cleanup=True):
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

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##################################################
# ensure that two jobs have two different names! #
##################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    jid = submit(local_buff_script, **kwargs)

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

def execBashCmdOnAllMachines(command, machines=MACHINES, mail=None, log=LOG_FILE):
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

    """
   https://docs.python.org/2/library/logging.html
   The logging module is intended to be thread-safe without any special work needing to be done by its clients.
   It achieves this though using threading locks; there is one lock to serialize access to the module’s shared data,
   and each handler also creates a lock to serialize access to its underlying I/O.

   If you are implementing asynchronous signal handlers using the signal module, you may not be able to use logging
   from within such handlers. This is because lock implementations in the threading module are not always re-entrant,
   and so cannot be invoked from such signal handlers.
        """

class condorThread(threading.Thread):
    def __init__(self, command=None, monitoringTimeStep=2, name=None, callBack=None, log=LOG_FILE, verbose=False):
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

    def logErr(self, msg):
        if self.verbose:
            logging.debug(msg)

    def run(self):
        self.logErr('run')
        if self.command:
            self.jid = submitWithBuffer(self.command, self.tName,
                                        mail=None,
                                        # group name of the job
                                        jobGroup='<group>',
                                        niceUser=False,
                                        requirements='(Memory > 1024)',
                                        priority=0,
                                        # maximum number of simultaneous jobs with the same group name
                                        maxSimultaneousJobsInGroup=100,
                                        log=self.log)
            while True:
                self.logErr('Monitoring')
                time.sleep(self.monitoringTimeStep)
                try:
                    (self.outFileName, self.errFileName) = getOutputWithBuffer(self.jid, self.tName, waitMaxTime=0, log=self.log)
                    self.res = (self.outFileName, self.errFileName)
                    self.logErr("thread %s has finished" % self.tName)
                    if self.callBack:
                        self.callBack(self.res)
                    break
                except NotYetFinished:
                    self.logErr('NotYetFinished')
                    continue
                except Exception as e:
                    self.logErr(str(type(e)) + ' ' + str(e))
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
    def example1_noBuff(nbJobs):
        listOfJNamesIds = []
        for jobName in range(nbJobs):
            # jid = job id (cluster ID of the new job)
            jid = submit("echo", arguments="Hello job %s!" % jobName)
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

    def example1_buff(nbJobs):
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

    def example1_thread(nbJobs):
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
    def example2_buff(nbJobs):
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

    def example2_thread(nbJobs):
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

    #######################################
    # Execute commands on remote machines #
    #######################################
    def example3():
        # echo 'toto
        execBashCmdOnAllMachines('echo "toto"')

        execBashCmdOnAllMachines('touch ' + REMOTE_BUFF_FOLDER + '/toto')
        # look inside the remote temporary folders to verify that temporary files have been removed
        execBashCmdOnAllMachines('ls ' + REMOTE_BUFF_FOLDER)

        # completely remove REMOTE_BUFF_FOLDER on all machines
        # If you have pbs on your jobs, the directory REMOTE_BUFF_FOLDER won't be clean automatically.
        # this script allows to remove all files in REMOTE_BUFF_FOLDER on all machines
        execBashCmdOnAllMachines('rm ' + REMOTE_BUFF_FOLDER + '/*')

    import timeit

    nbJobs = 20
    t_example1_noBuff = timeit.timeit("example1_noBuff(%s)" % nbJobs, setup="from __main__ import example1_noBuff", number=1)
    print >> sys.stderr, "example1_noBuff", t_example1_noBuff
    t_example1_buff = timeit.timeit("example1_buff(%s)" % nbJobs, setup="from __main__ import example1_buff", number=1)
    print >> sys.stderr, "example1_buff", t_example1_buff
    t_example1_thread = timeit.timeit("example1_thread(%s)" % nbJobs, setup="from __main__ import example1_thread", number=1)
    print >> sys.stderr, "t_example1_thread", t_example1_thread

    print >> sys.stderr, "example1_noBuff", t_example1_noBuff
    print >> sys.stderr, "example1_buff", t_example1_buff
    print >> sys.stderr, "t_example1_thread", t_example1_thread

    nbJobs = 20
    t_example2_buff = timeit.timeit("example2_buff(%s)" % nbJobs, setup="from __main__ import example2_buff", number=1)
    print >> sys.stderr, "t_example2_buff", t_example2_buff
    t_example2_thread = timeit.timeit("example2_thread(%s)" % nbJobs, setup="from __main__ import example2_thread", number=1)
    print >> sys.stderr, "t_example2_thread", t_example2_thread

    print >> sys.stderr, "t_example2_buff", t_example2_buff
    print >> sys.stderr, "t_example2_thread", t_example2_thread