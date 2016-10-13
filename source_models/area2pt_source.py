"""Read and OpenQuake area source model input file and convert
to a point source model, then write the point source file.

Jonathan Griffin
Geoscience Australia
October 2016
"""

import os, sys
import numpy
import hmtk

from openquake.hazardlib.source import area, point
#from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter, area_to_point_sources
from openquake.commonlib.sourcewriter import write_source_model
from openquake.commonlib.node import Node
from openquake.commonlib.sourcewriter import obj_to_node
from openquake.commonlib import nrml

def area2pt_source(area_source_file):
    """Calls OpenQuake parsers to read area source model
    from source_mode.xml type file, convert to point sources
    and write to a new nrml source model file.
    :params area_source_file:
        nrml format file of the area source
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=200.)
    parser = SourceModelParser(converter)
#    print [method for method in dir(parser) if callable(getattr(parser, method))]
    sources = parser.parse_sources(area_source_file) # will update at 2.1 to parse_group
    print sources

    name = 'test_point_model'
    #   source_model = Node("sourceModel", {"name": name}, nodes=nodes)
    nodes = []
    for source in sources:
        pt_sources = area_to_point_sources(source)
        nodes += (map(obj_to_node, pt_sources))
#    print nodes
    source_model = Node("sourceModel", {"name": name}, nodes=nodes)
    nrml_pt_file = area_source_file[:-4] + '_pts.xml'
    with open(nrml_pt_file, 'wb') as f:
        nrml.write([source_model], f, '%s')
#        for pt in pt_sources:
#nodes


if __name__ == "__main__":
    path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_inputs/NSHM'
    area_source_file = os.path.join(path, 'source_model_leonard_2008.xml')
    area2pt_source(area_source_file)
