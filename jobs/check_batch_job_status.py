#!/usr/bin/env python
"""
Script used for batch deaggregation jobs to check results were created
in each deagg procudction folder.  
Run script after the batch jobs script has been submitted and jobs have been
completed.  

To save information to disk use " python check_batch_job_status.py > tmp"
to save a temporary file to disk.  Then use "more tmp | grep 'FAIL'" to 
list failed jobs.  
"""

import os, sys
import re
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
