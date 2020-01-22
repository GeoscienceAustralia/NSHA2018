#!/usr/bin/env python
"""Wrapper script for running model and capturing provenance
"""

'''Edited by js1626 08/08/2019 to remove copying of source lt files'''


import os, sys
from os.path import join
from shutil import copy2
import getpass
import datetime
import argparse
import re
import csv
import subprocess

def find_replace(rep, infile_s, outfile_s):
    '''function to enter correct details in to job and parameter file '''
    
    infile = open(infile_s, "r")
    outfile = open(outfile_s, "wt")
    f_string = infile.read() 
    rep = dict((re.escape(k), v) for k, v in rep.iteritems())
    pattern = re.compile("|".join(rep.keys()))
    new_text = pattern.sub(lambda m: rep[re.escape(m.group(0))], f_string)
    outfile.write(new_text)
    
    infile.close()
    outfile.close()

# Set up input files from sites csv.

SA = "0.2"
SA_s = "SA02"
poe = "0.02"

with open("cities_021019.csv") as csvfile:
    csv_reader = csv.reader(csvfile, delimiter=',')
    job_file_list = []
    param_file_list = []
    
    for row in csv_reader:
        city = row[0]
        lat = row[1]
        lon = row[2]
        
        rep = {"<CITY>": city, "<SA>": SA, 
                "<LON>": lon, "<LAT>": lat, 
                "<POE>": poe, "<SA_s>": SA_s}

        ini_file = "job_deag_mll_TEMPLATE.ini"
        txt_file = "params_deag_mll_hi_mem_TEMPLATE.txt"
        ini_out = "job_deag_mll_" + str(city) + "_" + str(SA_s) + ".ini"
        txt_out = "params_deag_mll_hi_mem_" + str(city) +".txt"
        
        find_replace(rep, ini_file, ini_out)
        find_replace(rep, txt_file, txt_out)
        
        job_file_list.append(ini_out)
        param_file_list.append(txt_out)
        

# loop through the params file to set up directories and stuff?
output_dirs = []
user = getpass.getuser()
for i,param_file in enumerate(param_file_list):  

    params = {}
    try:
        f_in = open(param_file, 'r')
    except TypeError:
        #parser.print_help()
        raise Exception(('params file %s not found' % param_file))
    for line in f_in.readlines():
        row = line.split('=')
        params[row[0].rstrip()] = row[1].rstrip().lstrip()

    run_start_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    
    model_name = params['model_rel_path'].split('/')[-2]
    deag_name = params['deag_rel_path'].split('/')[-2]
    model_path = join(params['sandpit_path'], user, params['model_rel_path'])
    deag_path = join(params['sandpit_path'], user, params['deag_rel_path'])
    
    job_file = join(deag_path, params['job_file'])

    # Make output directory and copy input files
    model_output_base_city = params['model_output_base'][:-6]
    if not os.path.exists(model_output_base_city):
        os.mkdir(model_output_base_city)
    
    if not os.path.exists(params['model_output_base']):
        os.mkdir(params['model_output_base'])
             
    model_output_base_source_folder = join(params['model_output_base'], model_name)
    if not os.path.exists(model_output_base_source_folder):
        os.mkdir(model_output_base_source_folder)
        
    output_dir_name = '_'.join((params['job_type'], run_start_time, model_name, user))
    output_dir = join(model_output_base_source_folder, output_dir_name)
    os.mkdir(output_dir)
    
    # Copy files to output directory
    copy2(job_file, output_dir) #copy job file still required

    try:
        copy2(sites_file, output_dir)
    except IOError:
        msg = 'Sites file %s not found!' % sites_file
        raise IOError(msg)
    except NameError:
        print 'Warning: No sites file specified, is this intended?'

    # Build run_<model>.sh
    outlines = '#PBS -P w84\n'
    outlines += '#PBS -q normalbw\n' # for high-memory jobs
    outlines += '#PBS -l walltime=%s\n' % params['walltime']
    outlines += '#PBS -l ncpus=%s\n' % params['ncpus']
    outlines += '#PBS -l mem=%s\n' % params['mem']
    outlines += '#PBS -l wd\n'
    outlines += '#PBS -N oq512c512ht\n'
    outlines += '#PBS -l jobfs=%s\n' % params['jobfs']
    outlines += '#PBS -l other=hyperthread\n\n'
    
    #outlines += 'module load openquake/2.1.1\n'
    #outlines += 'module load openquake/2.4\n'
    #outlines += 'module load openquake/3.1\n' # used for NSHA18
    #outlines += 'module load openquake/3.3.1\n'
    outlines += 'module load openquake/3.6\n'
    outlines += 'oq-ini.all.sh\n'
    outlines += 'oq engine --run %s --exports csv >&  parjob.log\n' % params['job_file']
    outlines += 'oq-end.sh\n'

    run_script_name = 'run_%s.sh' % model_name
    run_script = join(output_dir, run_script_name)
    
    f_out = open(run_script, 'w')
    f_out.write(outlines)
    f_out.close()

# Change to output directory and submit job
# batch jobs need to go into ALL the folders and run the coresponding scripts
# there must be a way to loop through akl folder smartly.  
    output_dirs.append(output_dir)

for i,directory in enumerate(output_dirs):
    print directory
    os.chdir(directory)
    print(os.getcwd())
    cmd = 'qsub -v PBS_ARRAY_INDEX=%s %s' % (i+1,run_script_name)
    print cmd
    os.system(cmd)
    
# add section to clean up job files putting them in the right directories
