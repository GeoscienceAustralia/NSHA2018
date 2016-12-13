"""Read and OpenQuake area source model input file and convert
to a point source model, then write the point source file.

Jonathan Griffin
Geoscience Australia
October 2016
"""

import os, sys
import numpy
import hmtk
from subprocess import call

from openquake.hazardlib.source import area, point
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup
from openquake.commonlib.sourcewriter import write_source_model
from openquake.commonlib.node import Node
from openquake.commonlib.sourcewriter import obj_to_node
from openquake.commonlib import nrml

def area2pt_source(area_source_file, discretisation=200.):
    """Calls OpenQuake parsers to read area source model
    from source_mode.xml type file, convert to point sources
    and write to a new nrml source model file.
    :params area_source_file:
        nrml format file of the area source
    :params discretisation:
        Grid size (km) for the area source discretisation, 
        which defines the distance between resulting point
        sources.
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=discretisation)
    parser = SourceModelParser(converter)
    print [method for method in dir(parser)]# if callable(getattr(parser, method))]
    try:
        sources = parser.parse_sources(area_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(area_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    name = 'test_point_model'
    new_pt_sources = {}
    for source in sources:
        pt_sources = area_to_point_sources(source)
        for pt in pt_sources:
            pt.source_id = pt.source_id.replace(':','')
            pt.name = pt.name.replace(':','_')
            try:
                new_pt_sources[pt.tectonic_region_type].append(pt)
            except KeyError:
                new_pt_sources[pt.tectonic_region_type] = [pt]
           # print [method for method in dir(pt) if callable(getattr(pt, method))]
          #  print [attribute for attribute in dir(pt)]
    nrml_pt_file = area_source_file[:-4] + '_pts.xml'
    source_group_list = []
    id = 0
    for trt, sources in new_pt_sources.iteritems():
        source_group = SourceGroup(trt, sources = sources, id=id)
        id +=1
        source_group_list.append(source_group)
    write_source_model(nrml_pt_file, source_group_list,
                       name = 'Leonard2008')


if __name__ == "__main__":
    path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_inputs/NSHM'
    path = '/media/sf_openquake_shared_files/Australia/eq_hazard_tools/catalogue/'
    area_source_file = os.path.join(path, 'source_model_leonard_2008.xml')
    area2pt_source(area_source_file)
