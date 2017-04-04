import shapefile
from os import path
from numpy import array
from shapely.geometry import Point, Polygon
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

###############################################################################
# parse Leonard shp exported from OQ
###############################################################################

leoshp = path.join('shapefiles', 'NSHA13_regional_source_model.shp')

print 'Reading source shapefile...'
sf = shapefile.Reader(leoshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get shp data
code  = get_field_data(sf, 'name', 'str')
mmax  = get_field_data(sf, 'max_mag', 'float')
mmin  = get_field_data(sf, 'min_mag', 'float')
trt   = get_field_data(sf, 'trt', 'str')
dep_b = get_field_data(sf, 'hd1', 'float')

###############################################################################
# parse Leonard lookup csv to get completeness info
###############################################################################
l08_lookup = 'NSHA13_source_model.csv'

names = []
codes = []

# get codes and source names
lines = open(l08_lookup).readlines()[1:]

# zone numbers based on 2013 update
for line in lines:
    dat = line.strip().split(',')    
    if int(dat[0]) < 10:
        codes.append('Z00' + str(dat[0]))
        names.append('ZONE ' + str(dat[0]))
    elif int(dat[0]) >= 10 and int(dat[0]) < 100:
        codes.append('Z0' + str(dat[0]))
        names.append('ZONE ' + str(dat[0]))
    else:
        codes.append('Z' + str(dat[0]))
        names.append(dat[1])
    
###############################################################################
# get neotectonic domain number from centroid
###############################################################################
# load domains shp
dsf = shapefile.Reader(path.join('..','Domains','shapefiles','DOMAINS_NSHA18.shp'))
# get domains
neo_doms  = get_field_data(dsf, 'DOMAIN', 'float')
dom_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')

# get domain polygons
dom_shapes = dsf.shapes()
n_dom = []
n_mmax = []

# loop through L08 zones
for poly, code in zip(shapes, codes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    print clon, clat
    tmp_dom = 7
    tmp_mmax  = 8.0
    
    # loop through domains and find point in poly
    for neo_dom, mmax, dom_shape in zip(neo_doms, dom_mmax, dom_shapes):
        dom_poly = Polygon(dom_shape.points)
        
        # check if NSHA18 centroid in domains poly
        if point.within(dom_poly):
            tmp_dom = neo_dom
            tmp_mmax = mmax
            
        if code == 'Z019' or code == 'Z017':
            tmp_dom = 2
            tmp_mmax = 7.5

    n_dom.append(tmp_dom)
    n_mmax.append(tmp_mmax)

###############################################################################
# get TRT, depth, mcomp & ycomp form Leonard08
###############################################################################
# load domains shp
lsf = shapefile.Reader(path.join('..','Leonard2008','shapefiles','LEONARD08_NSHA18.shp'))

# get domains
ltrt  = get_field_data(lsf, 'TRT', 'str')
ldep  = get_field_data(lsf, 'DEP_BEST', 'float')
lycomp = get_field_data(lsf, 'YCOMP', 'str')
lmcomp = get_field_data(lsf, 'MCOMP', 'str')

# get domain polygons
l08_shapes = lsf.shapes()
trt = []
dep_b = []
ycomp = []
mcomp = []

# loop through L08 zones
for code, poly in zip(codes, shapes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    tmp_trt = -99
    tmp_dep = -99
    tmp_mc = -99
    tmp_yc = -99
    
    if code == 'Z036':
        tmp_trt = 'Non_cratonic'
        tmp_dep = 10
        tmp_mc = lmcomp[0]
        tmp_yc = lycomp[0]
        
    elif code == 'Z101' or code == 'Z102' or code == 'Z103':
        tmp_trt = 'Banda_Sea'
        tmp_dep = 10
        
        # use ISC-GEM completeness (Storchak et al., PEPI, 2015, doi: 10.1016/j.pepi.2014.06.009)
        tmp_yc = '1960;1918;1900'
        tmp_mc = '5.5;6.25;7.5'
        
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
    
    trt.append(tmp_trt)
    dep_b.append(tmp_dep)
    mcomp.append(tmp_mc)
    ycomp.append(tmp_yc)
    
dep_b = array(dep_b)             
   
###############################################################################
# write initial shapefile
###############################################################################

outshp = path.join('shapefiles', 'NSHA13_NSHA18.shp')

# set shapefile to write to
w = shapefile.Writer(shapefile.POLYGON)
w.field('SRC_NAME','C','50')
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
w.field('BVAL_BEST','F', 8, 5)
w.field('BVAL_LOWER','F', 8, 5)
w.field('BVAL_UPPER','F', 8, 5)
w.field('BVAL_FIX','F', 8, 2)
w.field('BVAL_FIX_S','F', 8, 2)
w.field('YCOMP','C','70')
w.field('MCOMP','C','50')
w.field('YMAX','F', 8, 0)
w.field('TRT','C','100')
w.field('DOMAIN','I', 2, 0)
w.field('CAT_FILE','C','50')

src_wt = 1.0
src_ty = 'area'
#dep_b = 10.
#dep_u = 5.
#dep_l = 15.
#min_mag = 4.8
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
#ycomp = '1980;1970;1965;1962;1958;1910;1880'
#mcomp = '3.0;3.5;4.0;4.5;5.0;6.0;6.4'
ymax  = 2016
#dom   = -99
cat   = 'GGcat-161025.csv'

# loop through original records
for i, shape in enumerate(shapes):

    # write new records
    if i >= 0:
        # # ignore banda sea sources
        if codes[i] == 'Z101' or codes[i] == 'Z102' or codes[i] == 'Z103':
            print 'Ignoring', codes[i]
    
        else:
            # set shape polygon
            w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
        
            dep_u = dep_b[i] - 0.5*dep_b[i]
            dep_l = dep_b[i] + 0.5*dep_b[i]
            w.record(names[i], codes[i], src_ty, str('%01d' % n_dom[i]), src_wt, dep_b[i], dep_u, dep_l, mmin[i], min_rmag, n_mmax[i], n_mmax[i]-0.2, n_mmax[i]+0.2, \
                     n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix, bval_fix_sig, ycomp[i], mcomp[i], ymax, trt[i], n_dom[i], cat)
        
# now save area shapefile
w.save(outshp)

# write projection file
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()
