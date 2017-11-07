"""Merge point sources built in parallel
"""

import os, sys
import time
import numpy as np
from NSHA2018.source_models.utils.pt2fault_distance import read_pt_source, \
    read_simplefault_source, combine_pt_sources
from glob import glob

source_model_name = 'National_Fault_Source_Model_2018_Collapsed_NSHA13'
area_source_model = '../zones/2012_mw_ge_4.0/NSHA13/input/collapsed/NSHA13_collapsed.xml'
#area_source_model = '../zones/2012_mw_ge_4.0/AUS6/input/collapsed/AUS6_collapsed.xml'
geom_pt_sources_filename =  area_source_model[:-4] + '_pts_geom_weighted.xml'

tmp_pt_source_filenames = glob(geom_pt_sources_filename.rstrip('.xml') + '_*.xml')
num_files = len(tmp_pt_source_filenames)/2
tmp_pt_source_filename_list = []
tmp_pt_source_list = []
# Now combine into one file
for j in range(0, num_files, 1):
    tmp_pt_filename = geom_filtered_pt_sources_sublist = geom_pt_sources_filename.rstrip('.xml') + \
        '_%03d.xml' % j 
    tmp_pt_source_filename_list.append(tmp_pt_filename)
for tmp_pt_source_file in tmp_pt_source_filename_list:
    print 'Reading %s' % tmp_pt_source_file
    tmp_pt_source = read_pt_source(tmp_pt_source_file)
    tmp_pt_source_list.append(tmp_pt_source)
merged_filename = geom_pt_sources_filename.rstrip('.xml') + '_merged_parallel.xml'
model_name = geom_pt_sources_filename.rstrip('.xml')
combine_pt_sources(tmp_pt_source_list, merged_filename, model_name, nrml_version = '04',
                   id_location_flag=None)


