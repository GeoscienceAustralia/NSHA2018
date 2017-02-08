"""Wrapper script for running model and capturing provenance
"""

import os, sys
from os.path import join
from subprocess import call
from shutil import copy2
import subprocess
import getpass
import datetime

# Parameters for building run_<model>.sh qsub job
# Normally you would only need to change these, although
# other parameters are hard coded at the bottom of the file
ncpus = 256
mem = '512GB'
walltime = '03:00:00'

# Paths to sandpit and input path, assuming you have
# a sandpit named after your NCI username
sandpit_path = '/short/w84/NSHA18/sandpit/'
model_rel_path = 'NSHA2018/source_models/smoothed_seismicity/Cuthbertson2016'
# All outputs will be under this directory
model_output_base = '/short/w84/NSHA18/PSHA_modelling/'

########################################################
# Nothing should change below here
####################################################### 

user = getpass.getuser()
run_start_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
model_name = model_rel_path.split('/')[-1]
model_path = join(sandpit_path, user, model_rel_path)
job_file = join(model_path, 'job.ini')

# Read job.ini file and find relevant input files
f_in = open(job_file, 'r')
for line in f_in.readlines():
    if line.startswith('source_model_logic_tree_file'):
        src_lt = line.split('=')[1].strip()
        src_lt_file = join(model_path, src_lt)
    if line.startswith('gsim_logic_tree_file'):
        gsim_lt = line.split('=')[1].strip()
        gsim_lt_file = join(model_path, gsim_lt)
    if line.startswith('export_dir'):
        export_dir = line.split('=')[1].strip()
f_in.close()
# Find source models from logic tree file
f_in = open(src_lt_file)
for line in f_in.readlines():
    if 'source_model' in line:
        source_model_name = line.replace('>','<').split('<')[2]
        print source_model_name
        source_model_file = join(model_path, source_model_name)        
f_in.close()    
# Make output directory and copy input files
output_dir_name = run_start_time + '_' + model_name + '_' + user
output_dir = join(model_output_base, output_dir_name)
os.mkdir(output_dir)
copy2(job_file, output_dir)
copy2(src_lt_file, output_dir)
copy2(gsim_lt_file, output_dir)
copy2(source_model_file, output_dir)

# Build run_<model>.sh
outlines = '#PBS -P w84\n'
outlines += '#PBS -q normal\n'
outlines += '#PBS -l walltime=%s\n' %walltime
outlines += '#PBS -l ncpus=%i\n' % ncpus
outlines += '#PBS -l mem=%s\n' % mem
outlines += '#PBS -l wd\n'
outlines += '#PBS -N oq512c512ht\n'
outlines += '#PBS -l jobfs=200GB\n'
outlines += '#PBS -l other=hyperthread\n\n'

outlines += 'module load openquake/2.1.1\n'
outlines += 'oq-ini.all.sh\n'
outlines += 'oq engine --run job.ini --exports csv >&  parjob.log\n'
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
