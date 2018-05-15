"""Utility function for parsing and converting area source models

Jonathan Griffin
Geoscience Australia April 2017
"""

import os, sys
import copy
import numpy as np
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

def weighted_pt_source(pt_sources, weights, name,
                       filename=None, nrml_version='04'):
    """Scales rates by weights for collapsing logic trees
    :param pt_sources:
        list of PointSource objects
    :param weights:
        dict contains weights for each tectonic
        region type, e.g. weights[trt] = 0.2
    :param filenam:
        path to output file, if provided will
        be written to nrml format as defined by
        nrml_version
    :param  nrml_version:
        version of nrml schema to use
    :returns weighted_pt_sources:
        list of PointSource objects with activity rates
        scaled by weights
    """
    weighted_point_sources = []
    for pt in pt_sources:
        new_pt = copy.deepcopy(pt) # Copy sources to avoid messing with original data
        mfd_type = type(pt.mfd).__name__
        trt = pt.tectonic_region_type
        weight = weights[trt]
        if mfd_type == 'TruncatedGRMFD':
            b_val = pt.mfd.b_val
            # rescale a value in log sapce
            a_val = np.log10(np.power(10, pt.mfd.a_val)*weight)
            new_pt.mfd.modify_set_ab(a_val, b_val)
        elif mfd_type == 'EvenlyDiscretizedMFD':
            mag_bins, rates = zip(*pt.mfd.get_annual_occurrence_rates())
            mag_bins = np.array(mag_bins)
            rates = np.array(rates)
            new_rates = rates*weight
            new_pt.mfd.modify_set_mfd(new_pt.mfd.min_mag, new_pt.mfd.bin_width,
                                      list(new_rates))
        else:
            msg = 'Weighting method for mfd type %s not yet defined' % mfd_type
            raise(msg)
        weighted_point_sources.append(new_pt)
    # Now write out
    if filename is not None:
        source_model_file = filename 
        print 'Writing to source model file %s' % source_model_file 
        if nrml_version == '04':
#            source_list = []
#            for source in weighted_point_sources:
#                source_list.append(source)
            nodes = list(map(obj_to_node, sorted(weighted_point_sources)))
            source_model = Node("sourceModel", {"name": name}, nodes=nodes)
            with open(source_model_file, 'wb') as f:
                nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
        elif nrml_version == '05':
            msg = 'Method not yet implemented for nrml version 0.5'
            raise(msg)
#            source_group_list = []
#            id = 0
#            for trt, sources in weighted_sources.iteritems():
#                source_group = SourceGroup(trt, sources = sources, id=id)
#                id +=1
#                source_group_list.append(source_group)
#            write_source_model(nrml_pt_file, source_group_list,
#                               name = name)
        else:
            print 'Warning: nrml version not specfied, xml not created'
    return weighted_point_sources

def area2pt_source(area_source_file, sources = None, investigation_time=50, 
                   rupture_mesh_spacing=10., width_of_mfd_bin=0.1,
                   area_source_discretisation=10.,
                   filename = None, nrml_version = '04',
                   name = None, return_faults = False,
                   exclude_ids = []):
    """Calls OpenQuake parsers to read area source model
    from source_mode.xml type file, convert to point sources
    and write to a new nrml source model file.
    :params area_source_file:
        nrml format file of the area source
    :params discretisation:
        Grid size (km) for the area source discretisation, 
        which defines the distance between resulting point
        sources.
    :params exclude_ids:
        Don't convert sources from given list of source_ids
        into points, preserve as area sources
    """
    if sources is None:
        sources = nrml2sourcelist(area_source_file, investigation_time=investigation_time, 
                                  rupture_mesh_spacing=rupture_mesh_spacing, 
                                  width_of_mfd_bin=width_of_mfd_bin,
                                  area_source_discretisation=area_source_discretisation)
    if name is None:
        name = '%s_points' % filename
    new_pt_sources = {}
    faults =[]
    excluded_area_sources = []
    for source in sources:
        if type(source).__name__ == 'ComplexFaultSource' or \
                 type(source).__name__ == 'SimpleFaultSource':
            faults.append(source)
        elif source.source_id in exclude_ids:
            excluded_area_sources.append(source)
        else:
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
    if return_faults:
        if len(exclude_ids) > 0:
            return source_group_list, faults, excluded_area_sources
        else:
            return source_group_list, faults
    else:
        if len(exclude_ids) > 0:
            return source_group_list, excluded_area_sources
        else:
            return source_group_list

