"""Converts an nrml file to a shapefile

Jonathan Griffin
Geoscience Australia
November 2016
"""

import os, sys
import numpy

from openquake.commonlib.shapefileparser import SourceModelParser, \
    ShapefileParser

def nrml2shp(nrml_source_file):
    """Converts an nrml file to a shapefile
    with attributes
    :params area_source_file:
        nrml file for an area based source model
    """
    shp_file_name = nrml_source_file[:-3] + 'shp'
    print shp_file_name
    parser = SourceModelParser()
    source_model = parser.read(nrml_source_file)
    shp_parser = ShapefileParser()
    shp_parser.write(shp_file_name, source_model)

if __name__ == "__main__":
    try:
        nrml_source_file = sys.argv[1]
    except:
        print 'Usage: python nrml2shp.py <source_model_nrml_file.xml>'
        sys.exit()
    nrml2shp(nrml_source_file)
