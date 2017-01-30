# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:20:55 2017

@author: u56903
"""

# gets the centroid (in lon/lat) from a list of polygon points
def get_shp_centroid(shp_points):
    from shapely.geometry import Polygon
    
    centroid = Polygon(shp_points).centroid.wkt
    clon, clat = centroid.strip('PIONT').replace('(',' ').replace(')',' ').split()
    
    return float(clon), float(clat)

# returns index of field name = field
def get_field_index(sf, field):
    fields = sf.fields[1:] # ignore first empty field
    for i, f in enumerate(fields):
        if f[0] == field:
            findex = i
            
    return findex

# returns a shapefile data field of name=field
def get_field_data(sf, field, datatype):
    '''
    sf = pyshp shapefile object
    field = data field name (str)    
    datatype = float, str
    
    returns:
    
    data = list of data in string or float fmt
    '''
    from misc_tools import checkfloat, checkint 
    
    # get index
    findex = get_field_index(sf, field)
    
    # get records
    recs = sf.records()
    
    # now loop thru recs and get data
    data = []
    for rec in recs:
        if datatype == 'str':
            data.append(rec[findex])
        elif datatype == 'float':
            data.append(checkfloat(rec[findex]))
        elif datatype == 'int':
            data.append(checkint(rec[findex]))
            
    return data
