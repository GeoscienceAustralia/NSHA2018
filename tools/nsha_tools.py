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
    
def get_shapely_centroid(poly):
    
    centroid = poly.centroid.wkt
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
    from numpy import array
    
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
            
    return array(data)


# gets preferred catalogue - uses NSHA cat if all vertices inside GG_cat_polygon
def get_preferred_catalogue(targetshpfile):
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data
    
    ###############################################################################
    # parse shapefile
    ###############################################################################
    
    catshpfile = '..//Other//ggcat_extent.shp'
    
    csf = shapefile.Reader(catshpfile)
    catshape = csf.shapes()[0] # only one shape
        
    tsf = shapefile.Reader(targetshpfile)
    targetshapes = tsf.shapes()
    
    ###############################################################################
    # find poly points within other polygons
    ###############################################################################
    
    # loop through target zones
    cat = []
    for poly in targetshapes:
        
        # set default catalogue
        tmpcat = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
        
        # now loop through all points in target shape
        for tlon, tlat in poly.points:
            point = Point(tlon, tlat)
            
            # check if point in catshape
            if point.within(Polygon(catshape.points)) == False:
                tmpcat = 'ISC-GEM_hmtk_declustered.csv'
                
        # now append catalogue
        cat.append(tmpcat)
        
    return cat

# converts datetime object to decmial years
# Slightly edited from: http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years    
def toYearFraction(date):
    from datetime import datetime as dt
    import time
    
    try:
        def sinceEpoch(date): # returns seconds since epoch
            return time.mktime(date.timetuple())
        s = sinceEpoch
    
        year = date.year
        startOfThisYear = dt(year=year, month=1, day=1)
        startOfNextYear = dt(year=year+1, month=1, day=1)
    
        yearElapsed = s(date) - s(startOfThisYear)
        yearDuration = s(startOfNextYear) - s(startOfThisYear)
        fraction = yearElapsed/yearDuration
    
        return date.year + fraction
    
    # for dates < 1900, work out manually to nearest month (should do something better, but fit for purpose!)
    except:
        return date.year + date.month/12.

# checks if string can be converted to float value.  If not returns nan        
def checkfloat(floatstr):
    try:
        return float(floatstr)
    except:
        from numpy import nan
        return nan

# checks if string can be converted to integer value.  If not returns nan       
def checkint(intstr):
    try:
        return int(intstr)
    except:
        from numpy import nan
        return nan   

def bval2beta(bval):
    from numpy import log
    return log(10**bval)

def beta2bval(beta):
    from numpy import log10, exp
    return log10(exp(beta))

