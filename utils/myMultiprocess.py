#!/usr/bin/python
# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Besides function available_cpu_count(), which is under CC by-SA 3.0 (see comment in function), all code is
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import os
import re
import subprocess

import sys
# For parallel computation
from multiprocessing import Queue, Process, Manager, Lock


# Matthieu's queue manager for parallel computations
def myPool(nbThreads, func, largs):

    # Librairies
    import multiprocessing
    import cStringIO

    # Sauvegarde des descripteurs
    backstdout = sys.stdout
    backstderr = sys.stderr

    # Semaphores
    sys.stdout = sys.stderr
    manager = multiprocessing.Manager()
    sys.stdout = backstdout
    queue = manager.Queue()

    def newfunc(i, args):
        queue.put((i, func(*args), sys.stdout.getvalue(), sys.stderr.getvalue()))

    def joinnext():
        (p, r, out, err) = queue.get()
        proc.pop(p).join()
        backstdout.write(out)
        backstderr.write(err)
        return (p, r)

    try:
        nrun = 0
        proc = {}
        for (i, x) in enumerate(largs):

            if nrun == nbThreads:
                yield joinnext()
                nrun -= 1

            # Lancement d'un nouveau processus
            proc[i] = multiprocessing.Process(target=newfunc, args=(i, x))
            sys.stdout = cStringIO.StringIO()
            sys.stderr = cStringIO.StringIO()
            proc[i].start()
            sys.stdout = backstdout
            sys.stderr = backstderr
            nrun += 1

        while len(proc) > 0:
            yield joinnext()

    finally:
        sys.stdout = backstdout
        sys.stderr = backstderr


# Wrapper used to print the % of the task progression
#######################################################
def wrapper(function):
    def f(input, kwargs, output, NbOfTasks, listOfPercentage, lock):
        for args in iter(input.get, 'STOP'):
            result = function(*args, **kwargs)
            output.put(result)
            progress = 100-100*(input.qsize()) / NbOfTasks
            lock.acquire()
            if progress in listOfPercentage:
                print >> sys.stderr, "%s" % progress + '%'
                listOfPercentage.remove(progress)
            lock.release()
    return f


# Multiprocess queue
#####################
# function is a reference to the function that will execute task's argsX.
# tasks = [(arg1Task1, arg2Task1), (arg2Task2, arg2Task2)]
# kwargs are the optional arguments of the function to be launched in parallel
def multiprocessTasks(function, tasks, **kwargs):
    #if kargs is not None:
    #    if 'verbose' in kargs:
    #        verbose = kargs['verbose']
    manager = Manager()
    # Lock is used to give priority to one process in order that it is not
    # disturbed by other processes.
    lock = Lock()
    # Manager is used to manage objects which are shared among processes.
    listOfPercentage = manager.list()
    for i in [j for j in range(0, 101, 5)[1:]]:
        listOfPercentage.append(i)

    #NUMBER_OF_PROCESSES = multiprocessing.cpu_count() * 2
    NUMBER_OF_PROCESSES = available_cpu_count()
    # TASKS = [ARGS-TAKS1, ARGS-TASK2, ...]
    # e.g. TASKS = [(c1, c2, g1[c1], g2[c2]) for (c1, c2) in itertools.product([c1 for c1 in g1], [c2 for c2 in g2])]
    TASKS = tasks
    # e.g. KWARGS = {'gapMax':gapMax, 'distanceMetric':distanceMetric, 'consistentSwDType':consistentSwDType, 'verbose':False}
    KWARGS = kwargs

    # Create queues
    task_queue = Queue()
    done_queue = Queue()

    # Submit tasks
    for task in TASKS:
        task_queue.put(task)

    # Start worker processes
    print >> sys.stderr, "synteny block extraction"
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=wrapper(function),
                args=(task_queue, KWARGS, done_queue, int(task_queue.qsize()),
                      listOfPercentage, lock)
                ).start()

    # Get and print results
    # 'Unordered results:'
    for i in range(len(TASKS)):
        for oneRes in done_queue.get():
            yield oneRes

    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')

    return


def available_cpu_count():
    # Function copied from Philipp Hagemeister
    # https://stackoverflow.com/questions/1006289
    """ Number of available virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling
    userspace-only program"""

    # cpuset
    # cpuset may restrict the number of *available* processors
    try:
        m = re.search(r'(?m)^Cpus_allowed:\s*(.*)$',
                      open('/proc/self/status').read())
        if m:
            res = bin(int(m.group(1).replace(',', ''), 16)).count('1')
            if res > 0:
                return res
    except IOError:
        pass

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # http://code.google.com/p/psutil/
    try:
        import psutil
        return psutil.NUM_CPUS
    except (ImportError, AttributeError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))

        if res > 0:
            return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])

        if res > 0:
            return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0:
            return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                  stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)

        if res > 0:
            return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')

        if res > 0:
            return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        res = 0
        for pd in pseudoDevices:
            if re.match(r'^cpuid@[0-9]+$', pd):
                res += 1

        if res > 0:
            return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0:
            return res
    except OSError:
        pass

    raise Exception('Can not determine number of CPUs on this system')
