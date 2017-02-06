# -*- coding: utf-8 -*-
"""
Created on Mon Feb 06 11:32:01 2017

@author: u56903
"""
import shapefile
from os import path
from numpy import array
from shapely.geometry import Point, Polygon
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid, get_shapely_centroid
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

###############################################################################
# parse Sinadinovski & McCue 2016 sources
###############################################################################
    
lines = open('sinadinovski_mccue.zones').readlines()

names = []
codes = []
polys = []
rawpolys = []
poly = []

for line in lines:
    if line.startswith('ZN'):
        codes.append(line.strip().split(':')[0].strip())
        names.append(line.strip().split(':')[1].strip())
        if len(poly) > 0:
            # close poly
            #poly.append(poly[0])
            # append to polys
            polys.append(Polygon(poly))
            rawpolys.append(poly)
        
        # reset poly
        poly = []
        
    # get ploy coords
    else:
        dat = line.strip().split(' ')
        poly.append([float(dat[1]), float(dat[0])])

# append last poly
#poly.append(poly[0])
polys.append(Polygon(poly))
rawpolys.append(poly)
#shapes = array(polygons)

###############################################################################
# get neotectonic domain number from centroid
###############################################################################
# load domains shp
dsf = shapefile.Reader(path.join('..','Domains','shapefiles','DOMAINS_NSHA18.shp'))

# get domains
neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')

# get domain polygons
dom_shapes = dsf.shapes()
dom = []
mmax = []

# loop through AUS6 zones
for code, poly in zip(codes, polys):
    # get centroid of leonard sources
    clon, clat = get_shapely_centroid(poly)
    point = Point(clon, clat)
    print clon, clat
    tmp_dom = -99
    tmp_mmax = -99
    
    # loop through domains and find point in poly    
    for neo_dom, neo_mx, dom_shape in zip(neo_doms, neo_mmax, dom_shapes):
        dom_poly = Polygon(dom_shape.points)
        
        # check if AUS6 centroid in domains poly
        if point.within(dom_poly):
            tmp_dom = neo_dom
            tmp_mmax = neo_mx
    
    dom.append(tmp_dom)
    mmax.append(tmp_mmax)
    
###############################################################################
# get TRT and depth form Leonard08
###############################################################################
# load domains shp
lsf = shapefile.Reader(path.join('..','Leonard2008','shapefiles','LEONARD08_NSHA18.shp'))

# get domains
ltrt  = get_field_data(lsf, 'TRT', 'str')
ldep  = get_field_data(lsf, 'DEP_BEST', 'float')

# get domain polygons
l08_shapes = lsf.shapes()
trt = []
dep_b = []

# loop through L08 zones
for poly in polys:
    # get centroid of leonard sources
    clon, clat = get_shapely_centroid(poly)
    point = Point(clon, clat)
    tmp_trt = -99
    tmp_dep = -99
    
    # loop through Leonard zones and find point in poly
    for zone_trt, zone_dep, l_shape in zip(ltrt, ldep, l08_shapes):
        l_poly = Polygon(l_shape.points)
        
        # check if leonard centroid in domains poly
        if point.within(l_poly):
            tmp_trt = zone_trt
            tmp_dep = zone_dep
    
    trt.append(tmp_trt)
    dep_b.append(tmp_dep)

dep_b = array(dep_b)    
###############################################################################
# write initial shapefile
###############################################################################

outshp = path.join('shapefiles','SIN_MCC_NSHA18_PUB.shp')

# set shapefile to write to
w = shapefile.Writer(shapefile.POLYGON)
w.field('SRC_NAME','C','100')
w.field('CODE','C','10')
#w.field('SRC_REGION','C','100')
#w.field('SRC_REG_WT','F', 8, 3)
w.field('SRC_TYPE','C','10')
w.field('SRC_WEIGHT','F', 8, 2)
w.field('DEP_BEST','F', 8, 1)
w.field('DEP_UPPER','F', 8, 1)
w.field('DEP_LOWER','F', 8, 1)
w.field('MIN_MAG','F', 8, 2)
w.field('MIN_RMAG','F', 8, 2)
w.field('MMAX_BEST','F', 8, 2)
w.field('MMAX_LOWER','F', 8, 2)
w.field('MMAX_UPPER','F', 8, 2)
w.field('N0_BEST','F', 8, 5)
w.field('N0_LOWER','F', 8, 5)
w.field('N0_UPPER','F', 8, 5)
w.field('BVAL_BEST','F', 8, 5)
w.field('BVAL_LOWER','F', 8, 5)
w.field('BVAL_UPPER','F', 8, 5)
w.field('BVAL_FIX','F', 8, 2)
w.field('BVAL_FIX_S','F', 8, 2)
w.field('YCOMP','C','70')
w.field('MCOMP','C','30')
w.field('YMAX','F', 8, 0)
w.field('TRT','C','100')
w.field('DOMAIN','F', 2, 0)
w.field('CAT_FILE','C','50')

src_wt = 1.0
src_ty = 'area'
#dep_b = 10.
dep_u = dep_b - 0.5*dep_b
dep_l = dep_b + 0.5*dep_b
min_mag = 4.8
min_rmag = 2.5
#mmax[i]
#mmax_l = mmax[i]-0.2
#mmax_u = mmax[i]+0.2
n0 = -99
n0_l = -99
n0_u = -99
bval = -99
bval_l = -99
bval_u = -99
bval_fix = -99
bval_fix_sig = -99
ycomp = '1880;1910;1958;1962;1965;1970;1980'
mcomp = '6.4;6.0;5.0;4.5;4.0;3.5;3.0'
ycomp = '1980;1970;1965;1962;1958;1910;1880'
mcomp = '3.0;3.5;4.0;4.5;5.0;6.0;6.4'
ymax  = 2016
#trt   = 'TBD'
#dom   = -99
cat   = 'GGcat-161025.csv'

# loop through original records
for i, rawpoly in enumerate(rawpolys):

    # set shape polygon
    w.line(parts=[rawpoly], shapeType=shapefile.POLYGON)
        
    # write new records
    if i >= 0:
        w.record(names[i], codes[i], src_ty, src_wt, dep_b[i], dep_u[i], dep_l[i], min_mag, min_rmag, mmax[i], mmax[i]-0.2, mmax[i]+0.2, \
                 n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix, bval_fix_sig, ycomp, mcomp, ymax, trt[i], dom[i], cat)
        
# now save area shapefile
w.save(outshp)

# write projection file
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()       