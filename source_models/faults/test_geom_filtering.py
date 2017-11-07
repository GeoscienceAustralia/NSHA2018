"""Run test datasets for geometrically filtered approach
"""

import os, sys
import numpy as np
from NSHA2018.source_models.utils.area_sources import nrml2sourcelist, area2pt_source
from NSHA2018.source_models.utils.pt2fault_distance import read_simplefault_source, \
    pt2fault_distance

# Basic parameters 
filepath = './Adelaide'
fsm = os.path.join(filepath, 'Adelaide_faults.xml')
area_source_model = os.path.join(filepath, 'source_model_adelaide.xml')
bin_width = 0.1
investigation_time = 50
fault_mesh_spacing = 2 #2 Fault source mesh
rupture_mesh_spacing = 2 #10 # Area source mesh
area_source_discretisation = 10 #20

# Read in the area source model
print 'Reading area source model %s' % area_source_model
area_pt_filename = area_source_model[:-4] + '_pts.xml'
area_sources = nrml2sourcelist(area_source_model, 
                               investigation_time=investigation_time, 
                               rupture_mesh_spacing=rupture_mesh_spacing, 
                               width_of_mfd_bin=bin_width,
                               area_source_discretisation=area_source_discretisation)
# Convert area sources to point sources for filtering
print 'Converting to point sources'
name = area_source_model.split('/')[-1][:-4] + '_pts' 
point_sources = area2pt_source(area_source_model, sources=area_sources,
                               filename=area_pt_filename,
                               name = name)
pt_source_list = []
for source_group in point_sources:
     for source in source_group:
          pt_source_list.append(source)

print 'Applying geometrical filtering'
name = area_source_model.split('/')[-1][:-4] + '_pts_geom_filter' 
fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)
revised_point_sources = area_source_model[:-4] + '_pts_geom_filtered.xml'
pt2fault_distance(pt_source_list, fault_sources, min_distance = 21.0,
                  filename = revised_point_sources,
                  buffer_distance = 100.,
                  name = name)

#print 'Combining revised pt sources with fault source model into one file'
