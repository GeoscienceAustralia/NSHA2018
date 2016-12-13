"""Determines whether a point is within a polygon, and depending
on the answer performs an action

Jonathan Griffin
Geoscience Australia
November 2016
"""

import os, sys
import shapely
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter#, SourceGroup
from openquake.commonlib.sourcewriter import write_source_model

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

def read_area_source(area_source_file, discretisation=200.):
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
    for source in sources:
        print source.polygon
    return sources

def edit_pt_in_polygon(pt_source, polygon, attribute):
    """Determines which polygon a point source is in, and adjusts the 
    relevant attribute based on the polygon
    """
    pass

if __name__ == "__main__":
    pt_source_file = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_runs/smoothed_seismicity_adaptive/source_model_smoothed_adaptive_K4_0.1.xml'
    area_source_file = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_runs/NSHM_background/source_model_NSHM_background.xml'
    attribute = 'tectonicRegion'

#    pt_sources = read_pt_source(pt_source_file)
    area_sources = read_area_source(area_source_file)
    edit_pt_in_polygon(pt_sources, area_sources, attribute)
