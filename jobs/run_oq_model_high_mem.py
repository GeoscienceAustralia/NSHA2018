#!/usr/bin/env python2
"""Wrapper script for running model and capturing provenance
"""

import os, sys
from os.path import join
from subprocess import call
from shutil import copy2
import subprocess
import getpass
import datetime
import argparse

parser = argparse.ArgumentParser(description='Build NCI run scripts for OpenQuake models' \
                                     'capturing appropriate provenance data')
parser.add_argument('-param_file', type=str, default=None,
                    help='Name of file containing model run parameters')
args = parser.parse_args()

params = {}
try:
    f_in = open(args.param_file, 'r')
except TypeError:
    parser.print_help()
    raise Exception(('params file %s not found' % args.param_file))
for line in f_in.readlines():
    row = line.split('=')
    params[row[0].rstrip()] = row[1].rstrip().lstrip()

user = getpass.getuser()
run_start_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
model_name = params['model_rel_path'].split('/')[-1]
model_path = join(params['sandpit_path'], user, params['model_rel_path'])
job_file = join(model_path, params['job_file'])

# Make output directory and copy input files
if not os.path.exists(params['model_output_base']):
    os.mkdir(params['model_output_base'])
model_output_base_source_folder = join(params['model_output_base'], model_name)
if not os.path.exists(model_output_base_source_folder):
    os.mkdir(model_output_base_source_folder)
output_dir_name = '_'.join((params['job_type'], run_start_time, model_name, user))
output_dir = join(model_output_base_source_folder, output_dir_name)
print 'Output Path:', output_dir
os.mkdir(output_dir)

# make full job path
full_job_dir = join(output_dir, params['job_file'])
print full_job_dir

# We want to ensure we are using the same ground motion models
# for all model runs
gsim_lt_file = join(params['sandpit_path'], user, params['shared_rel_path'],
                    params['gsim_lt_filename'])

# Read job.ini file and find relevant input files
f_in = open(job_file, 'r')
for line in f_in.readlines():
    if line.startswith('source_model_logic_tree_file'):
        src_lt = line.split('=')[1].strip()
        src_lt_file = join(model_path, src_lt)
    if line.startswith('gsim_logic_tree_file'):
        gsim_lt = line.split('=')[1].strip()
        if gsim_lt != params['gsim_lt_filename']:
            msg = 'Ground motion logic tree filename %s as specified in run_oq_model.py' \
                'different to %s specified in %s' \
                % (params['gsim_lt_filename'], gsim_lt, job_file.split('/')[-1])
            raise ValueError(msg)
    if line.startswith('sites_csv'):
        sites_filename = line.split('=')[1].strip()
        sites_file = join(params['sandpit_path'], user, params['shared_rel_path'],
                          sites_filename)
    if line.startswith('export_dir'):
        export_dir = line.split('=')[1].strip()
f_in.close()

# Copy files to output directory
copy2(job_file, output_dir)
copy2(src_lt_file, output_dir)
copy2(gsim_lt_file, output_dir)
#copy2('job_stats.ini', output_dir)
try:
    copy2(sites_file, output_dir)
except IOError:
    msg = 'Sites file %s not found!' % sites_file
    raise IOError(msg)
except NameError:
    print 'Warning: No sites file specified, is this intended?'

# Find source models from logic tree file
# and copy to output dir
f_in = open(src_lt_file)
for line in f_in.readlines():
    if '.xml' in line:
        source_model_name = line.replace('>','<').split('<')[2]
        source_model_file = join(model_path, source_model_name)        
        copy2(source_model_file, output_dir)
f_in.close()    

# Build run_<model>.sh
outlines = '#PBS -P w84\n'
outlines += '#PBS -q normal\n' # for high-memory jobs
outlines += '#PBS -l storage=normal/w84\n'
outlines += '#PBS -l walltime=%s\n' % params['walltime']
outlines += '#PBS -l ncpus=%s\n' % params['ncpus']
outlines += '#PBS -l mem=%s\n' % params['mem']
outlines += '#PBS -l wd\n'
outlines += '#PBS -N oq512c512ht\n'
outlines += '#PBS -l jobfs=%s\n' % params['jobfs']
outlines += '#PBS -l other=hyperthread\n\n'

#outlines += 'module load openquake/2.1.1\n'
#outlines += 'module load openquake/2.4\n'
#outlines += 'module load openquake/3.1\n' # used in the NSHA18
#outlines += 'module load openquake/3.6\n'
outlines += 'module load openquake/3.7.1\n'
outlines += 'oq-ini.all.sh\n'
# zip inputs - oq zip /path/to/your/job.ini job.zip
#outlines += 'oq zip %s job.zip\n' % full_job_dir
outlines += 'oq info %s >& info.log\n' % params['job_file']
outlines += 'oq engine --run %s --exports csv >&  parjob.log\n' % params['job_file']
#outlines += 'oq engine --lhc >&  lhc.log\n'
#outlines += 'oq engine --run job_stats.ini --hc 1 --exports csv >&  jobstats.log\n'
outlines += 'oq-end.sh'

run_script_name = 'run_%s.sh' % model_name
run_script = join(output_dir, run_script_name)
f_out = open(run_script, 'w')
f_out.write(outlines)
f_out.close()

# Change to output directory and submit job
os.chdir(output_dir)
cmd = 'qsub %s' % run_script_name
print cmd
os.system(cmd)
