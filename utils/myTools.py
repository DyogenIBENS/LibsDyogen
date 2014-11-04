# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import os
import sys
import itertools
import time
import string
import warnings

import enum

from functools import wraps

import myFile

null = open('/dev/null', 'w')

debug = null

class Namespace: pass

def applyFunctions(fun, data):
    for (f, x) in itertools.izip(fun, data):
        yield f(x)

def funcFilter(fun):
    return lambda data: (f(x) for (f, x) in itertools.izip(fun, data))

def __delitem__(self, key):
    dict.__delitem__(self, self[key])
    dict.__delitem__(self, key)

def __len__(self):
    """Returns the number of connections"""
    return dict.__len__(self) // 2

# Print a correctly spaced table, all columns are aligned
# Input table format:
# table =[[l0c0, l0c1, l0c2],
#         [l1c0, l1c1, l1c2]]
# With 'lxcy' a string that corresponds to line 'x' and column 'y'
def printTable(table, output):
    max_lens = []
    for i in range(len(table[0])):
        max_lens.append(max([len(str(r[i])) for r in table]))
    res = "\n".join(
        ["".join([string.ljust(str(e), l + 2)
                  for e, l in zip(r, max_lens)])
         for r in table])
    print >> output, res
    return res


# decorator that adds a switchable verbose mode to a function
# FIXME, if the decorated function is called with more arguments that in the
# function definition, there is no error raised.
def verbose(functionToExcecute):
    @wraps(functionToExcecute) # to avoid changing the name of the function
    def modifiedFunction(*args, **kargs):
        if 'verbose' in kargs:
            if kargs['verbose'] == True:
                res = functionToExcecute(*args, **kargs)
                # **kargs still contains verbose
            else:
                sys.stderr = open(os.devnull, 'w')
                res = functionToExcecute(*args, **kargs)
                # **kargs still contains verbose
                sys.stderr = sys.__stderr__
        else:
            warnings.warn("function %s has no option verbose although it uses a verbose decorator" % functionToExcecute.__name__, category=SyntaxWarning, stacklevel=2)
            res = functionToExcecute(*args, **kargs)
        return res
    return  modifiedFunction

# decorator for functions that requires a minimal python version >= 2.7 for instance
# version is a tuple, for instance if the function requires python version at least 2.7, version = (2,7)
def minimalPythonVersion(version):
    def decorator(functionToExcecute):
        def modifiedFunction(*args, **kargs):
            if sys.version_info < version:
                raise Exception("Function %s needs at least python %s.%s" % (functionToExcecute.__name__,version[0],version[1]))
            else:
                return functionToExcecute(*args, **kargs)
        return modifiedFunction
    return decorator

# decorator that computes the execution time
def tictac(functionToExcecute):
    @wraps(functionToExcecute) # to avoid changing the name of the function
    def modifiedFunction(*args,**kargs):
        tic = time.time()
        res = functionToExcecute(*args,**kargs)
        tac = time.time()
        deltaTicTac = tac - tic
        print >> sys.stderr, "Function \"%s\" was executed in %s seconds" % (functionToExcecute.__name__, deltaTicTac)
        return res
    return modifiedFunction

# decorator that warns the user that the function is deprecated
def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emmitted
    when the function is used."""
    import warnings
    def newFunc(*args, **kwargs):
        warnings.warn("Call to deprecated function %s." % func.__name__, category=DeprecationWarning, stacklevel=2)
        return func(*args, **kwargs)
    newFunc.__name__ = func.__name__
    newFunc.__doc__ = func.__doc__
    newFunc.__dict__.update(func.__dict__)
    return newFunc


# record results of a function for each parameter value #
class memoize:
    """Decorator that caches a value returned by the function each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """
    def __init__(self, func):
        self.func = func
        self.nbcall = 0
        self.cache = {}

    def __repr__(self):
        return "[%s: %d values cached, %d calls, %.2fx speedup]" %\
            (self.func.__name__, len(self.cache),
             self.nbcall,
             self.nbcall/float(len(self.cache)) if len(self.cache) > 0 else 0)

    def __call__(self, *args, **kwargs):
        self.nbcall += 1
        try:
            return self.cache[args]
        except KeyError:
            self.cache[args] = self.func(*args)
            value = self.func(*args)
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            print >> sys.stderr, "Warning: %s is not cacheable (from %s/%s)" % (args, self, self.func)
            return self.func(*args)

    def reinit_stats(self):
        self.nbcall = 0
        self.cache = {}

    def __doc__(self):
        """Return the function's docstring."""
        return self.func.__doc__

class memoizeMethod(object):
    def __init__(self, function):
        self.function = function
        self.nbcall = 0
        self.cache = {}

    def __get__(self, instance, cls=None):
        self.instance = instance
        return self

    def __call__(self, *args):
        self.nbcall += 1
        if args in self.cache:
            return self.cache[args]
        else:
            res = self.cache[args] = self.function(self.instance, *args)
            return res

    def reinit_stats(self):
        self.nbcall = 0
        self.cache = {}


# FIXME, best memoizeMethod
# cf http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/
#class memoizeMethod2(object):
#    """cache the return value of a method
#
#
#    This class is meant to be used as a decorator of methods. The return value
#    from a given method invocation will be cached on the instance whose method
#    was invoked. All arguments passed to a method decorated with memoize must
#    be hashable.
#
#    If a memoized method is invoked directly on its class the result will not
#    be cached. Instead the method will be invoked like a static method:
#    class Obj(object):
#        @memoize
#        def add_to(self, arg):
#            return self + arg
#    Obj.add_to(1) # not enough arguments
#    Obj.add_to(1, 2) # returns 3, result is not cached
#    """

# Check is an excecutable is accessible
# This may be usefull to check if a plugged external programm has been added to
# the PATH environment variable.
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


# iterator of adjacent components of a list
class myIterator:
    # sliding couple (x,y)
    @staticmethod
    def slidingTuple(lst):
        if len(lst) > 0:
            x = lst[0]
            for i in xrange(1, len(lst)):
                y = lst[i]
                yield (x,y)
                x = y

# liste of partitions of size k in range(n)
@memoize
def partitions(n, k):
    if n == 1:
        if k == 1:
            return [ [[0]] ]
    if (n >= k) and (k >= 1):
        all = []
        for x in partitions(n-1, k-1):
            all.append( x + [[n-1]] )
        for x in partitions(n-1, k):
            for i in xrange(k):
                all.append( [y if i != j else y + [n-1] for (j,y) in enumerate(x)] )
        return all
    else:
        return []

# management of a parallel execution on a range of values
def getRange(s):
    if myFile.hasAccess(s):
        f = myFile.openFile(s, "r")
        lst = []
        for l in f:
            lst.extend( [int(x) for x in l.replace('\n', '').split()] )
        f.close()
        return lst
    else:
        (start,_,end) = s.partition(':')
        return range(int(start), int(end)+1)


# hashable dict class, useful to use it as a key
class hashabledict(dict):
    def __hash__(self):
        return hash(tuple(sorted(self.items())))

# hashable list class
class hashablelist(list):
    def __hash__(self):
        return hash(tuple(self))

# This class allows to group a list of elements.
# From an initial list of elements, links are added between these elements.
# The class gathers elements that are linked.
class myCombinator:

    def __init__(self, ini = []):
        self.grp = list(ini)
        self.dic = {}
        for i in xrange(len(self.grp)):
            self.grp[i] = list(set(self.grp[i]))
            for x in self.grp[i]:
                self.dic[x] = i

    # define a link between all elements of obj
    # update sets already built
    def addLink(self, obj):

        if len(obj) == 0:
            return []

        obj = set(obj)
        grp = self.grp
        dic = self.dic

        # elements of obj already present in the combinator
        d = set( dic[x] for x in obj if x in dic )

        if len(d) == 0:
            # None, the obj is added just like it is
            i = len(grp)
            grp.append(list(set(obj)))
            for x in obj:
                dic[x] = i
            return grp
        else:
            i = d.pop()
            grpiextend = grp[i].extend
            for x in d:
                grpiextend(grp[x])
                for y in grp[x]:
                    dic[y] = i
                #FIXME not del grp[x] ?
                # see reduce and "empty sets"
                grp[x] = []
            dd = [x for x in obj if x not in dic]
            for x in dd:
                dic[x] = i
            grpiextend(dd)
            return grp[i]


    # return an iterator over the data
    # empty sets are thus removed
    def __iter__(self):
        for g in self.grp:
            if len(g) > 0:
                yield g

    # remove empty sets
    def reduce(self):
        self.__init__(self)

    # reset combinator
    def reset(self):
        self.__init__()

# add more options to a specific module
__moduleoptions = []
def addModuleOptions(namespace, options):
    for (name,typ,val) in options:
        __moduleoptions.append( (namespace+":"+name,typ,val) )


# ask a list of file in arguments
class FileList:
    def __init__(self, value):
        self.minNbFiles = value

    def __repr__(self):
        return '<FileList(%d)>' % self.minNbFiles


# Parse arguments on the command line
#  1. requested arguments (name, builder)
#  2. options in the form of -opt=val (name, builder, default_value)
# If an error occurs, user's command line is printed as well as a short description of the bug and a brief manual of the script (info).
def checkArgs(args, options, info, showArgs=True):

    options = options + __moduleoptions
    # print error informations if wrong arguments
    def error_usage(reason):
        print >> sys.stderr, "- ERROR -", reason
        print >> sys.stderr, " Usage :", sys.argv[0]
        for (i,t) in enumerate(args):
            print >> sys.stderr, "\t", "%d:" % (i+1), t[0], t[1]
        for t in options:
            if isinstance(t[1], enum.Enum):
                print >> sys.stderr, "\t", "  -%s %s (%s)" % (t[0], t[1]._keys, t[2])
            elif t[1] == bool:
                print >> sys.stderr, "\t", "+/-%s (%s)" % (t[0],t[2])
            else:
                print >> sys.stderr, "\t", "  -%s %s (%s)" % t
        if info != "":
            print >> sys.stderr, "\n", info
        sys.exit(1)

    def putValue(typ, val, v):
        # instantiate the value depending on the type
        if typ == bool:
            # Type booleen
            res = {"false": False, "true":True}[v.lower()]
        elif typ == file:
            # Type 'fichier': test of presence
            v = os.path.expanduser(v)
            if not myFile.hasAccess(v):
                error_usage("File '%s' innaccessible" % v)
            else:
                res = v
        elif isinstance(typ, enum.Enum):
            try:
                res = getattr(typ, v)
            except AttributeError:
                error_usage("'%s' is not among %s" % (v,typ._keys))
        else:
            # otherwise the builder is used
            res = typ(v)
            if isinstance(val, list) and (res not in val):
                # non authorised parameter value
                error_usage("'%s' is not among %s" % (res,myFile.myTSV.printLine(val, '/')))
        return res

    valOpt = {}
    valArg = {}
    opt = {}
    for (name,typ,val) in options:
        opt[name] = (typ,val)
        valOpt[name] = val[0] if isinstance(val, list) else getattr(typ, val) if isinstance(typ, enum.Enum) else val

    # arguments are scanned, counted and values are extracted
    for tt in sys.argv[1:]:

        t = tt.replace('^', ' ')

        # an optional argument
        if t[0] in '-+':

            # non bool parameter
            try:
                i = t.index('=')
                s = t[1:i]

                # the parameter name must be known
                if not s in valOpt:
                    error_usage("Option '%s' unknown" % s)

                valOpt[s] = putValue(opt[s][0], opt[s][1], t[i+1:])

            # if '=' is not found, it is a bool type
            except ValueError:
                s = t[1:]
                # unexpected parameter name
                if s not in valOpt:

                    # predefined values
                    if s.startswith("psyco"):
                        if t[0] == '+':
                            try:
                                import psyco
                                psyco.full()
                            except ImportError:
                                print >> sys.stderr, "Unable to load psyco !"
                    elif s == "bz2":
                        if t[0] == '+':
                            import bz2
                            sys.stdout = bz2.BZ2File("/dev/stdout", "w")
                    elif s == "gz":
                        if t[0] == '+':
                            import gzip
                            sys.stdout = gzip.GzipFile("/dev/stdout", "w")
                    elif (s == "lzma") or (s == "xz"):
                        if t[0] == '+':
                            import lzma
                            sys.stdout = lzma.LZMAFile("/dev/stdout", "w")
                    elif s == "debug":
                        if t[0] == '+':
                            global debug
                            debug = sys.stderr
                    else:
                        error_usage("Option '%s' unknown" % s)
                elif opt[s][0] != bool:
                    error_usage("Use -%s=value" % s)
                else:
                    # Here, False is assigned
                    valOpt[s] = (t[0] == '+')
        else:
            if len(valArg) < len(args):
                (s,typ) = args[len(valArg)]
                if isinstance(typ, FileList):
                    valArg[s] = list()
                    assert len(valArg) == len(args)
                    valArg[s].append(putValue(file, None, t))
                else:
                    valArg[s] = putValue(typ, None, t)
            elif isinstance(args[-1][1], FileList):
                valArg[args[-1][0]].append(putValue(file, None, t))
            else:
                error_usage("Too many arguments on '%s'" % t)

    if isinstance(args[-1][1], FileList):
        if args[-1][0] not in valArg:
            valArg[args[-1][0]] = []
        if len(valArg[args[-1][0]]) < args[-1][1].minNbFiles:
            error_usage("Not enough files for '%s'" % args[-1][0])

    # there is less than the minimal number of arguments
    if len(valArg) < len(args):
        error_usage("Not enough arguments")

    valArg.update(valOpt)
    if showArgs:
        print >> sys.stderr, "Arguments:", valArg
    return valArg
