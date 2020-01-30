#!/usr/bin/env python

import os, sys
from os.path import join
from glob import glob

f = open("Job_list_test.txt", 'r')
for line in f:
    #print(line)
    print("**********")
    line = re.sub('\n', '', line).strip() # remove \n characters
    city = line.split("/")[6]
    full_line = join(line, "results")
    #print(full_line)
    if os.path.exists(full_line):
        print("PASS: results folder exists for %s" % city)
    else:
        print("FAIL: No results folder created for %s" % city)
    rlz_files = glob("%s/rlz*" %full_line)
    if len(rlz_files) == 3:
        print("PASS: 3 results files created for %s" % city)
    elif len(rlz_files) < 3:
        print("FAIL: less than 3 output files for %s" % city)
f.close()
