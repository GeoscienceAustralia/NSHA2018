"""Call geometrical filter in parallel to speed up
"""

import os, sys
import time
from time import localtime, strftime, gmtime
import string
import numpy as np
import pypar
from NSHA2018.source_models.utils.pt2fault_distance import read_pt_source, \
    read_simplefault_source, pt2fault_distance, combine_pt_sources

# Set up paralell
proc = pypar.size()                # Number of processors as specified by mpirun                     
myid = pypar.rank()                # Id of of this process (myid in [0, proc-1])                     
node = pypar.get_processor_name()  # Host name on which current process is running                   
print 'I am proc %d of %d on node %s' % (myid, proc, node)
#nruns = 320 # currently hard coded - need to improve this                                            
t0 = pypar.time()

fault_mesh_spacing = 2 #2 Fault source mesh                     
rupture_mesh_spacing = 2 #10 # Area source mesh                                                         
area_source_discretisation = 15 #20 
source_model_name = 'National_Fault_Source_Model_2018_Collapsed_DIMAUS_2018'
#area_source_model = '../zones/2018_mw/NSHA13/input/collapsed/NSHA13_collapsed.xml'
#area_source_model = '../zones/2018_mw/NSHA13/input/collapsed/NSHA13_collapsed.xml'
area_source_model = '../zones/2018_mw/DIMAUS/input/collapsed/DIMAUS_collapsed.xml'
geom_pt_sources_filename =  area_source_model[:-4] + '_pts_geom_weighted.xml'
geom_pt_sources = read_pt_source(geom_pt_sources_filename)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

# Split sources
list_length = len(geom_pt_sources) / (proc*10)
print list_length
if (len(geom_pt_sources) % proc) > 0:
    list_length +=1
pt_list = list(chunks(geom_pt_sources, list_length))
#print pt_list
print len(pt_list)
#sys.exit()
fsm =  os.path.join(source_model_name, source_model_name + '_geom_filtered.xml')
fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)

#pt_filename_list = []
for i in range(0, len(pt_list), 1):
    if i % proc == myid:
        run = "%03d" % i
        # Apply geometrical filtering                                                                        
        print 'Applying geometrical filtering for run %s' % run
        geom_filtered_pt_sources = geom_pt_sources_filename.rstrip('.xml') + \
            '_' + run + '.xml'
        try:
            pt2fault_distance(pt_list[i], fault_sources, min_distance=5.0,
                              filename=geom_filtered_pt_sources,
                              buffer_distance = 50.,
                              name=source_model_name)
        except IndexError:
            print 'List index %i out of range' % i


pypar.barrier()
if myid == 0:
    tmp_pt_source_filename_list = []
    tmp_pt_source_list = []
    # Now combine into one file
    for j in range(0, len(pt_list), 1):
        tmp_pt_filename = geom_filtered_pt_sources_sublist = geom_pt_sources_filename.rstrip('.xml') + \
            '_%03d.xml' % j 
        tmp_pt_source_filename_list.append(tmp_pt_filename)
    for tmp_pt_source_file in tmp_pt_source_filename_list:
        tmp_pt_source = read_pt_source(tmp_pt_source_file)
        tmp_pt_source_list.append(tmp_pt_source)
    merged_filename = geom_pt_sources_filename.rstrip('.xml') + '_merged_parallel.xml'
    model_name = geom_pt_sources_filename.rstrip('.xml')
    combine_pt_sources(tmp_pt_source_list, merged_filename, model_name, nrml_version = '04',
                       id_location_flag=None)
#
#if myid == 0:
    ss = int(pypar.time() - t0)
    h = ss / 3600
    m = (ss % 3600) / 60
    s = (ss % 3600) % 60
    print "--------------------------------------------------------"
    print 'P0: Total time (%i seconds): %s:%s:%s (hh:mm:ss)' % (ss,
                                                                string.zfill(h, 2),
                                                                string.zfill(m, 2),
                                                                string.zfill(s,2))
    print "--------------------------------------------------------"
pypar.finalize()
