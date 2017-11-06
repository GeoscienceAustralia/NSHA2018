"""Code to test how the pt sources are combined
"""

import os, sys
from NSHA2018.source_models.logic_trees import logic_tree
from NSHA2018.source_models.smoothed_seismicity.combine_ss_models import combine_ss_models
from NSHA2018.source_models.utils.pt2fault_distance import read_pt_source, combine_pt_sources
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.mfd.evenly_discretized import EvenlyDiscretizedMFD

infiles = ['Australia_Adaptive_K4_b1.198.xml', 'Australia_Adaptive_K4_b1.352.xml',
           'Australia_Adaptive_K4_b1.043.xml']
domains_shp = '../zones/2012_mw_ge_4.0/NSHA13_Background/shapefiles/NSHA13_BACKGROUND_NSHA18_M\
FD.shp'

# Create test datasets
create_test_data =False
if create_test_data == True:
    for pt_file in infiles:
        pt_sources = read_pt_source(pt_file)
        pt_sources = pt_sources[:50]
        name = 'test'
        nodes = list(map(obj_to_node, sorted(pt_sources)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        outfile = pt_file.rstrip('.xml') + '_testdata.xml'
        with open(outfile, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
if create_test_data == False:
    pass
lt  = logic_tree.LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
filedict_bestb = {'Non_cratonic': 'Australia_Adaptive_K4_b1.198_testdata.xml',
                  'Cratonic': 'Australia_Adaptive_K4_b1.198_testdata.xml',
                  'Extended': 'Australia_Adaptive_K4_b1.198_testdata.xml'}
filedict_upperb = {'Non_cratonic': 'Australia_Adaptive_K4_b1.043_testdata.xml',
                   'Cratonic': 'Australia_Adaptive_K4_b1.043_testdata.xml',
                   'Extended': 'Australia_Adaptive_K4_b1.043_testdata.xml'}
filedict_lowerb = {'Non_cratonic': 'Australia_Adaptive_K4_b1.352_testdata.xml',
                   'Cratonic': 'Australia_Adaptive_K4_b1.352_testdata.xml',
                   'Extended': 'Australia_Adaptive_K4_b1.352_testdata.xml'}

combine_ss_models(filedict_bestb, domains_shp, lt, 'test_bestb.xml', nrml_version = '04', weight=0.5)
combine_ss_models(filedict_upperb, domains_shp, lt, 'test_upperb.xml', nrml_version = '04', weight=0.3)
combine_ss_models(filedict_lowerb, domains_shp, lt, 'test_lowerb.xml', nrml_version = '04', weight=0.2)

# combine all pt source models                                                                 
point_source_list = ['test_bestb.xml', 'test_upperb.xml', 'test_upperb.xml']
filename = 'test_combine_merge.xml'
filepath = filename
name = filename.rstrip('.xml')

# read list of files                                                                           
pt_source_model_list =[]
for point_source_model in point_source_list:
    print 'Reading %s' % point_source_model
    pt_model = read_pt_source(point_source_model)
    pt_source_model_list.append(pt_model)
combine_pt_sources(pt_source_model_list, filepath, name , nrml_version='04',
                   id_location_flag = 'location')
merge_pts = read_pt_source(filepath)
print 'Final files has %i points' % (len(merge_pts))
