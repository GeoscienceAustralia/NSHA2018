"""Reads an input point source model and generates ruptures according 
to the MFD and nodal plane distirbution. Nearest distance to the fault is calculated and used to reduce the Mmax of the point source model. 
"""

import os, sys
import numpy as np
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter
from openquake.commonlib.sourcewriter import write_source_model
from openquake.commonlib.node import Node
from openquake.commonlib.sourcewriter import obj_to_node
from openquake.commonlib import nrml

def read_pt_source(pt_source_file):
    """Read nrml source model into pt source objects
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=200.)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(pt_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(pt_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    return sources

def read_simplefault_source(simplefault_source_file):
    """Read nrml source model into simpmle fault objects
    """
    converter = SourceConverter(50, 2, width_of_mfd_bin=0.1,
                                area_source_discretization=200.)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(simplefault_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(simplefault_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    return sources


if __name__ == "__main__":
    path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_inputs/NSHM'
    path = '/media/sf_openquake_shared_files/Australia/eq_hazard_tools/catalogue/'
    pt_source_file = os.path.join(path, 'source_model_leonard_2008_pts.xml')
    fault_path = '/home/openquake/GEM/eq_hazard_tools/shapefile2nrml/Adelaide_test/'
    simplefault_source_file = os.path.join(fault_path, 'Adelaide_faults.xml')
    pt_sources = read_pt_source(pt_source_file)
    fault_sources = read_simplefault_source(simplefault_source_file)
