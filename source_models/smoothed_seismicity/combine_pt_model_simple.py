"""This code combines the smoothed seismicity models applying different b-values and mmax's in different tectonic regions
"""

import os, sys
import numpy as np
import ogr
import shapefile
from shapely.geometry import Point, Polygon

from NSHA2018.source_models.logic_trees import logic_tree
from NSHA2018.source_models.utils.pt2fault_distance import read_pt_source, combine_pt_sources
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.mfd.evenly_discretized import EvenlyDiscretizedMFD



if __name__ == "__main__":
    outfile_bestb = 'source_model_Australia_Adaptive_K3_merged_bestb.xml'
    outfile_upperb = 'source_model_Australia_Adaptive_K3_merged_upperb.xml'
    outfile_lowerb = 'source_model_Australia_Adaptive_K3_merged_lowerb.xml'
    point_source_list = [outfile_bestb, outfile_upperb, outfile_lowerb]
    filename = 'source_model_Australia_Adaptive_K3_merged_inc_b_mmax_uncert_v1.xml'
    name = filename.rstrip('.xml')

    # read list of files
#    pt_source_model_list =[]
#    for point_source_model in point_source_list:
#        print 'Reading %s' % point_source_model
#        pt_model = read_pt_source(point_source_model)
##        pt_source_model_list.append(pt_model)
##    combine_pt_sources(pt_source_model_list, filename, name , nrml_version='04',
##                       id_location_flag = 'id')

    outfile_bestb = 'smoothed_frankel_50_3_mmin_3.0_merged_bestb.xml'
    outfile_upperb = 'smoothed_frankel_50_3_mmin_3.0_merged_upperb.xml'
    outfile_lowerb = 'smoothed_frankel_50_3_mmin_3.0_merged_lowerb.xml'
    point_source_list = [outfile_bestb, outfile_upperb, outfile_lowerb]
    filename = 'source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert_v1.xml'
    name = filename.rstrip('.xml')

    # read list of files
    pt_source_model_list =[]
    for point_source_model in point_source_list:
        print 'Reading %s' % point_source_model
        pt_model = read_pt_source(point_source_model)
        pt_source_model_list.append(pt_model)
    combine_pt_sources(pt_source_model_list, filename, name , nrml_version='04',
                       id_location_flag = 'id')

