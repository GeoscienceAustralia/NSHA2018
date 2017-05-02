"""Utility function for parsing and converting area source models

Jonathan Griffin
Geoscience Australia April 2017
"""

import os, sys

from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup
from openquake.hazardlib.sourcewriter import write_source_model, obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml

def nrml2sourcelist(area_source_file, investigation_time=50, 
                    rupture_mesh_spacing=10., width_of_mfd_bin=0.1,
                    area_source_discretisation=10.):
    """Parser nrml file containing area sources and read into
    a list of source objects
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=area_source_discretisation)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(area_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(area_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    return sources

def area2pt_source(area_source_file, sources = None, investigation_time=50, 
                   rupture_mesh_spacing=10., width_of_mfd_bin=0.1,
                   area_source_discretisation=10.,
                   filename = None, nrml_version = '04',
                   name = None):
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
    if sources is None:
        sources = nrml2sourcelist(area_source_file, investigation_time=investigation_time, 
                                  rupture_mesh_spacing=rupture_mesh_spacing, 
                                  width_of_mfd_bin=width_of_mfd_bin,
                                  area_source_discretisation=area_source_discretisation)
    if name is None:
        name = '%s_points' % filename
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
    if filename is not None:
        if nrml_version == '04':
            source_list = []
            for trt, sources in new_pt_sources.iteritems():
                for source in sources:
                    source_list.append(source)
            nodes = list(map(obj_to_node, sorted(source_list)))
            source_model = Node("sourceModel", {"name": name}, nodes=nodes)
            with open(nrml_pt_file, 'wb') as f:
                nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
        # This will write version 0.5
        elif nrml_version == '05':
            write_source_model(nrml_pt_file, source_group_list,
                               name = filename)
        else:
            print 'Warning: nrml version not specfied, xml not created'
    return source_group_list
