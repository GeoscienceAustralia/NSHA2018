"""Utilities functions for smoothed seismicity models
"""

import ogr
import shapefile
from shapely.geometry import Point, Polygon
import numpy as npy

def params_from_shp(shapefile):
    """Get parameters from shapefile attribute table
    """
    print 'Getting completeness and b-values from %s' % shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile, 0)
    dsf = data_source.GetLayer()
    param_list = []
    for feature in dsf:
        mcomp = feature.GetField('MCOMP').split(';')
        ycomp = feature.GetField('YCOMP').split(';')
        completeness_table = []
        for i in range(len(mcomp)):
            completeness_table.append([float(ycomp[i]),float(mcomp[i])])
        #completeness_table = np.array(completeness_table)
        params = {'BVAL_BEST': feature.GetField('BVAL_BEST'),
                  'BVAL_LOWER': feature.GetField('BVAL_LOWER'),
                  'BVAL_UPPER': feature.GetField('BVAL_UPPER'),
                  'COMPLETENESS': completeness_table}
        if not any(d == params for d in param_list):
            param_list.append(params)
    for params in param_list:
        params['COMPLETENESS'] = np.array(params['COMPLETENESS'])
    return param_list
