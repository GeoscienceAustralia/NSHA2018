"""Script to build smoothed seismicity model that only uses the additive method (i.e.
has full weight of 1 for the model)
"""

import os, sys
import numpy as np
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_simplefault, \
    shapefile_2_simplefault_CE, shapefile_2_simplefault_MM, shp2nrml
from NSHA2018.source_models.logic_trees import logic_tree
from NSHA2018.source_models.utils.utils import largest_remainder
from NSHA2018.source_models.utils.area_sources import nrml2sourcelist, \
    area2pt_source, weighted_pt_source
from NSHA2018.source_models.utils.pt2fault_distance import read_simplefault_source, \
    pt2fault_distance, write_combined_faults_points, combine_pt_sources, read_pt_source
#from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from subprocess import call

from openquake.hazardlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup



def build_model(combined_output_dir, pt_source_model, fsm,
                source_model_name, fault_mesh_spacing):
    """Combine the fault and pt source model"""

 #   additive_pt_sources_filename =  pt_source_model[:-4] + '_pts.xml'
    print 'reading fault sources'
    additive_pt_sources = read_pt_source(pt_source_model)
    # Apply additive approach                                                                                     
    if not os.path.exists(combined_output_dir):
        os.makedirs(combined_output_dir)
    print 'Writing full additive model'
    outfile =  os.path.join(combined_output_dir, source_model_name)
    fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)
    write_combined_faults_points(additive_pt_sources, fault_sources,
                                 outfile, source_model_name, nrml_version = '04')

if __name__ == "__main__":

    fsm = '../faults/National_Fault_Source_Model_2018_Collapsed_NSHA13/National_Fault_Source_Model_2018_Collapsed_NSHA13_all_methods_collapsed_inc_cluster.xml'
    fault_mesh_spacing = 2 #2 Fault source mesh
    combined_output_dir = 'GA_adaptive_smoothing_collapsed_faults_K5'
    ss_source_model = 'GA_adaptive_smoothing_collapsed/source_model_Australia_Adaptive_K5_merged_inc_b_mmax_uncert.xml'
    source_model_name = 'GA_adaptive_smoothing_collapsed_faults_additive'
    print 'Building Adaptive SS with faults'
    build_model(combined_output_dir, ss_source_model, fsm,
                source_model_name, fault_mesh_spacing)

    combined_output_dir = 'GA_fixed_smoothing_collapsed'
    pt_source_model = 'GA_fixed_smoothing_collapsed/source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert.xml'
    source_model_name = 'GA_fixed_smoothing_collapsed_faults_additive'
    print 'Building Fixed SS with faults'
    build_model(combined_output_dir, pt_source_model, fsm,
                source_model_name, fault_mesh_spacing)

