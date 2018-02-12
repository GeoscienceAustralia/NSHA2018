import shapefile
from os import path
from shapely.geometry import Point, Polygon
from numpy import array, zeros_like, where, median, std
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid, get_preferred_catalogue
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

###############################################################################
# parse Domains shp
###############################################################################

domshp = 'Domains_Sep2011_edit.shp'

print 'Reading source shapefile...'
sf = shapefile.Reader(domshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get src name
src_name = get_field_data(sf, 'CODE', 'str')
#src_domain = get_field_data(sf, 'Id', 'float')

###############################################################################
# parse Domains lookup csv
###############################################################################

domcsv = 'Domains_Sep2011_lookup.csv'

dom = []
mmax = []
codes = []
name = []
mcomp = []
ycomp = []

lines = open(domcsv).readlines()[1:]
for line in lines:
    dat = line.strip().split(',')
    dom.append(int(round(float(dat[0]))))
    codes.append(dat[1])
    '''
    name.append(dat[2])
    mmax.append(float(dat[3]))
    mcomp.append(dat[4])
    ycomp.append(dat[5])
    #trt.append(dat[6])
    '''

# reset Gawler Craton to Flinders due to b-value similarities
dom[4] = 2

# reset W Tasmania to Flinders due to b-value similarities (Delamerian Orogen)
dom[1] = 2


###############################################################################
# get preferred catalogues 
###############################################################################

# get preferred catalogues for each zone
prefCat = get_preferred_catalogue(domshp)
    
###############################################################################
# load 2018 completeness models
###############################################################################

# load domains shp
compshp = path.join('..','Other','Mcomp_NSHA18_smoothed.shp')
mcsf = shapefile.Reader(compshp)

# get completeness data
mc_ycomp = get_field_data(mcsf, 'YCOMP', 'str')
mc_mcomp = get_field_data(mcsf, 'MCOMP', 'str')

# get completeness polygons
mc_shapes = mcsf.shapes()

# set empty completeness values
ycomp = []
mcomp = []
min_rmag = []

# loop through Mcomp zones
for code, poly in zip(codes, shapes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    print clon, clat
    
    # loop through domains and find point in poly    
    matchidx = -99
    mccompFound = False
    for i in range(0, len(mc_shapes)):
        dom_poly = Polygon(mc_shapes[i].points)
        
        # check if Domains centroid in completeness poly
        if point.within(dom_poly): 
            ycomp.append(mc_ycomp[i])
            mcomp.append(mc_mcomp[i])
            mccompFound = True
    
    # if no Mcomp model assigned, use conservative model
    if mccompFound == False:
        ycomp.append('1980;1964;1900')
        mcomp.append('3.5;5.0;6.0')
        
    # set rmin range
    min_rmag.append(max([3.0, float(mcomp[-1].split(';')[0])]))
    
# use manual modification
min_rmag = [3.0, 3.0, 3.0, 3.5, 3.0, 3.0, 3.2, 3.5, 3.3, 3.3, 3.0, 3.3, 3.5, 3.5, 3.5, 3.3]

###############################################################################
# get neotectonic domains number and Mmax from zone centroid
###############################################################################
# get path to reference shapefile
#shapepath = open('..//reference_shp.txt').read() # this is used for other sources!
shapepath = 'Domains_NSHA18.shp'

#print '\nNOTE: Getting Domains info for original magnitudes\n'
#shapepath = open('..//reference_shp_mx.txt').read()

# load domains shp
dsf = shapefile.Reader(shapepath)

# get domains
neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
neo_bval = get_field_data(dsf, 'BVAL_BEST', 'float')
neo_bval_l = get_field_data(dsf, 'BVAL_LOWER', 'float')
neo_trt  = get_field_data(dsf, 'TRT', 'str')
neo_dep  = get_field_data(dsf, 'DEP_BEST', 'float')

# get bval sigma
bval_sig = neo_bval_l - neo_bval

# get domain polygons
dom_shapes = dsf.shapes()
#dom = []
mmax = []
trt = []
dep_b = []
bval_fix = []
bval_sig_fix = []

# loop through Domains zones
for code, poly in zip(codes, shapes):
    # get centroid of leonard sources
    clon, clat = get_shp_centroid(poly.points)
    point = Point(clon, clat)
    print clon, clat
        
    # loop through domains and find point in poly
    matchidx = -99
    for i in range(0, len(dom_shapes)):
        dom_poly = Polygon(dom_shapes[i].points)
        
        # check if Domain centroid in domains poly
        if point.within(dom_poly):
            matchidx = i
    
    if code == 'OSGB' or code == 'EAPM' or code == 'WARM':
        matchidx = -1
    #elif code == 'GAWL':
    #    matchidx = 0
            
    # set dummy values
    if matchidx == -99:
        dom.append(-99)
        mmax.append(-99)
        trt.append(-99)
        dep_b.append(-99)
        bval_fix.append(-99)
        bval_sig_fix.append(-99)
    
    # fill real values
    else:        
        #dom.append(neo_doms[matchidx])
        mmax.append(neo_mmax[matchidx])
        trt.append(neo_trt[matchidx])
        dep_b.append(neo_dep[matchidx])
        #ycomp.append(neo_ycomp[matchidx])
        #mcomp.append(neo_mcomp[matchidx])
        #bval_fix.append(neo_bval[matchidx])
        #bval_sig_fix.append(bval_sig[matchidx])
        print '\nNOTE: Setting b-value params to -99\n'
        print code, matchidx
        bval_fix.append(-99)
        bval_sig_fix.append(-99)
        
    

dep_b = array(dep_b)

###############################################################################
# load Rajabi SHMax vectors 
###############################################################################

shmaxshp = path.join('..','Other','SHMax_Rajabi_2016.shp')

print 'Reading SHmax shapefile...'
sf = shapefile.Reader(shmaxshp)
    
# get src name
shmax_lat = get_field_data(sf, 'LAT', 'float')
shmax_lon = get_field_data(sf, 'LON', 'float')
shmax     = get_field_data(sf, 'SHMAX', 'float')

###############################################################################
# get preferred strike
###############################################################################
shmax_pref = []
shmax_sig  = []

for code, poly in zip(codes, shapes):
    # get shmax points in polygon
    shm_in = []
    
    # now loop through earthquakes in cat
    for shmlo, shmla, shm in zip(shmax_lon, shmax_lat, shmax):
        
        # check if pt in poly and compile mag and years
        pt = Point(shmlo, shmla)
        if pt.within(Polygon(poly.points)):
            shm_in.append(shm)
    
    if len(shm_in) > 0: 
        shmax_pref.append(median(array(shm_in)))
        shmax_sig.append(std(array(shm_in)))
        print 'Getting SHmax for', code
    
    # if no points in polygons, get nearest neighbour
    else:
        print 'Getting nearest neighbour...'
        min_dist = 9999.
        for shmlo, shmla, shm in zip(shmax_lon, shmax_lat, shmax):
            pt = Point(shmlo, shmla)
            pt_dist = pt.distance(Polygon(poly.points))
            if pt_dist < min_dist:
                min_dist = pt_dist
                shm_near = shm
        
        shmax_pref.append(shm_near) # set nearest neighbour
        shmax_sig.append(15.) # set std manually
                 
###############################################################################
# write initial shapefile
###############################################################################

outshp = 'Domains_NSHA18.shp'

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
w.field('SHMAX','F', 6, 2)
w.field('SHMAX_SIG','F', 6, 2)
w.field('YMAX','F', 8, 3)
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
#min_rmag = 3.0
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
ymax  = 2017

# loop through original records
for i, shape in enumerate(shapes):

    # set shape polygon
    w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
        
    # write new records
    if i >= 0:
        w.record(src_name[i], codes[i], src_ty, dom[i], src_wt, dep_b[i], dep_u[i], dep_l[i], min_mag, min_rmag[i], mmax[i], mmax[i]-0.2, mmax[i]+0.2, \
                 n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix[i], bval_sig_fix[i], ycomp[i], mcomp[i], shmax_pref[i], shmax_sig[i], ymax, trt[i], dom[i], prefCat[i])
        
# now save area shapefile
w.save(outshp)

# write projection file
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()
