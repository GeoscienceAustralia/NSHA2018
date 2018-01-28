import shapefile
from os import path
from numpy import array, zeros_like, where
from shapely.geometry import Point, Polygon
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid, get_preferred_catalogue
except:
    print 'Add PYTHONPATH to NSHA18 root directory'


###############################################################################
# parse Mcomp shp
###############################################################################

#arupshp = 'ARUP_source_model.shp'
mcshp = 'completeness_zones.shp'

# get preferred catalogues
prefCat = get_preferred_catalogue(mcshp)

print 'Reading source shapefile...'
sf = shapefile.Reader(mcshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get src name
src_name = get_field_data(sf, 'NAME', 'str')
codes = src_name

###############################################################################
# parse ARUP lookup csv
###############################################################################
## input file modified to classify Tasmainia as non-cratonic from super domains
#arupcsv = 'ARUP_source_model.edit.csv'
#
#name = []
#codes = []
#neo_domains = []
#
#lines = open(arupcsv).readlines()[1:]
#for line in lines:
#    dat = line.strip().split(',')
#    name.append(dat[1])
#    codes.append('ZN'+dat[3]) # use "sub_zone" instead
#    neo_domains.append(dat[2])
#    
###############################################################################
# get neotectonic superdomains number and Mmax from zone centroid
###############################################################################
# get path to reference shapefile
shapepath = open('..//reference_shp.txt').read()

#print '\nNOTE: Getting Domains info for original magnitudes\n'
#shapepath = open('..//reference_shp_mx.txt').read()

# load domains shp
dsf = shapefile.Reader(shapepath)

# get domains
neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
neo_bval = get_field_data(dsf, 'BVAL_BEST', 'float')
neo_bval_l = get_field_data(dsf, 'BVAL_LOWER', 'float') # lower curve, higher b-value
neo_trt  = get_field_data(dsf, 'TRT', 'str')
neo_dep  = get_field_data(dsf, 'DEP_BEST', 'float')
neo_ycomp = get_field_data(dsf, 'YCOMP', 'str')
neo_mcomp = get_field_data(dsf, 'MCOMP', 'str')

# get bval sigma
bval_sig = neo_bval_l - neo_bval

# get domain polygons
dom_shapes = dsf.shapes()
dom = []
mmax = []
trt = []
nclass = []
dep_b = []
ycomp = []
mcomp = []
bval_fix = []
bval_sig_fix = []

# loop through ARUP zones
for code, poly in zip(codes, shapes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    print clon, clat
    
    # loop through domains and find point in poly    
    matchidx = -99
    for i in range(0, len(dom_shapes)):
        dom_poly = Polygon(dom_shapes[i].points)
        
        # check if ARUP centroid in domains poly
        if point.within(dom_poly): 
            matchidx = i
    
    if code.startswith('SA'):
        matchidx = 5
        print 'Fixing index: ', code
    elif code.startswith('ZN1b'):
        matchidx = 0
        print 'Fixing index: ', code
            
    # set dummy values
    if matchidx == -99:
        dom.append(-99)
        mmax.append(-99)
        trt.append(-99)
        dep_b.append(-99)
        ycomp.append(-99)
        mcomp.append(-99)
        bval_fix.append(-99)
        bval_sig_fix.append(-99)
        nclass.append(-99)
    # fill real values
    else:
        dom.append(neo_doms[matchidx])
        mmax.append(neo_mmax[matchidx])
        trt.append(neo_trt[matchidx])
        dep_b.append(neo_dep[matchidx])
        ycomp.append(neo_ycomp[matchidx])
        mcomp.append(neo_mcomp[matchidx])
        bval_fix.append(neo_bval[matchidx])
        bval_sig_fix.append(bval_sig[matchidx])
        #nclass.append(-99)

dep_b = array(dep_b)
    
"""    
###############################################################################
# get TRT, depth and completeness info from Leonard08
###############################################################################
# load Leonard08 shp
lsf = shapefile.Reader(path.join('..','Leonard2008','shapefiles','LEONARD08_NSHA18.shp'))

# get leonard08 data
ltrt  = get_field_data(lsf, 'TRT', 'str')
ldep  = get_field_data(lsf, 'DEP_BEST', 'float')
lycomp = get_field_data(lsf, 'YCOMP', 'str')
lmcomp = get_field_data(lsf, 'MCOMP', 'str')

# get leonard08 polygons
l08_shapes = lsf.shapes()

# set variables to be filled

# loop thru ARUP shapes
for code, poly in zip(codes, shapes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    tmp_trt = -99
    tmp_dep = -99
    tmp_mc = -99
    tmp_yc = -99
    
    # set Mmax values for zones outside of Leonard08
    if code == 'ZN7b':
        tmp_trt = 'Non_cratonic'
        tmp_dep = 10
        tmp_mc = lmcomp[0]
        tmp_yc = lycomp[0]
        
    elif code == 'MSAB':
        tmp_trt = 'Non_cratonic'
        tmp_dep = 10
        tmp_mc = lmcomp[0]
        tmp_yc = lycomp[0]
        
    # loop through Leonard zones and find point in poly
    else:
        for zone_trt, zone_dep, yc, mc, l_shape \
            in zip(ltrt, ldep, lycomp, lmcomp, l08_shapes):
            l_poly = Polygon(l_shape.points)
            
            # check if leonard centroid in domains poly
            if point.within(l_poly):
                tmp_trt = zone_trt
                tmp_dep = zone_dep
                tmp_mc = mc
                tmp_yc = yc
    
    # append new data
    trt.append(tmp_trt)
    dep_b.append(tmp_dep)
    mcomp.append(tmp_mc)
    ycomp.append(tmp_yc)

dep_b = array(dep_b)             
"""
###############################################################################
# write initial shapefile
###############################################################################

outshp = 'Mcomp_NSHA18.shp'

# set shapefile to write to
w = shapefile.Writer(shapefile.POLYGON)
w.field('SRC_NAME','C','100')
w.field('CODE','C','10')
#w.field('SRC_REGION','C','100')
#w.field('SRC_REG_WT','F', 8, 3)
w.field('SRC_TYPE','C','10')
w.field('CLASS','C','10')
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
w.field('BVAL_BEST','F', 8, 3)
w.field('BVAL_LOWER','F', 8, 3)
w.field('BVAL_UPPER','F', 8, 3)
w.field('BVAL_FIX','F', 8, 3)
w.field('BVAL_FIX_S','F', 8, 3)
w.field('YCOMP','C','70')
w.field('MCOMP','C','50')
w.field('YMAX','F', 8, 0)
w.field('TRT','C','100')
w.field('DOMAIN','F', 2, 0)
w.field('CAT_FILE','C','50')

src_wt = 1.0
src_ty = 'area'
#dep_b = 10.

dep_u = 0.5 * array(dep_b)

idx = where(dep_b > 7.)[0]
dep_l = zeros_like(dep_b)
dep_l[idx] = 1.5 * array(dep_b[idx])
idx = where(dep_b < 7.)[0]
dep_l[idx] = 2 * array(dep_b[idx])

min_mag = 4.5
min_rmag = 3.0
#mmax[i]
#mmax_l = mmax[i]-0.2
#mmax_u = mmax[i]+0.2
n0 = -99
n0_l = -99
n0_u = -99
bval = -99
bval_l = -99
bval_u = -99
#bval_fix = -99
#bval_fix_sig = -99
'''
ycomp = '1880;1910;1958;1962;1965;1970;1980'
mcomp = '6.4;6.0;5.0;4.5;4.0;3.5;3.0'
ycomp = '1980;1970;1965;1962;1958;1910;1880'
mcomp = '3.0;3.5;4.0;4.5;5.0;6.0;6.4'
'''
ymax  = 2011
#trt   = 'TBD'
#dom   = -99

# loop through original records
for i, shape in enumerate(shapes):

    # set shape polygon
    w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
        
    # write new records
    if i >= 0:
        w.record(codes[i], codes[i], src_ty, dom[i], src_wt, dep_b[i], dep_u[i], dep_l[i], min_mag, min_rmag, mmax[i], mmax[i]-0.2, mmax[i]+0.2, \
                 n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix[i], bval_sig_fix[i], ycomp[i], mcomp[i], ymax, trt[i], dom[i], prefCat[i])
        
# now save area shapefile
w.save(outshp)

# write projection file
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()
