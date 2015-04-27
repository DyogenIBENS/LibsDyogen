#!/bin/bash
# cythonise the utils/extractDiags.pyx file into utils/extractSbs.c
cython utils/extractSbs.pyx
# C compilation of utils/extractSbs.c into utils/extractSbs.so that can be used as a python package
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o utils/extractSbs.so utils/extractSbs.c