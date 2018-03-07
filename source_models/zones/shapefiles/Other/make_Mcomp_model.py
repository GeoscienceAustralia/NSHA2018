import shapefile
from shapely.geometry import Polygon
from numpy import ones_like, array
try:
    from tools.nsha_tools import get_field_data
    from tools.source_shapefile_builder import get_preferred_catalogue, \
                                               get_completeness_model, get_aus_shmax_vectors, \
                                               get_rate_adjust_factor, build_source_shape, \
                                               get_ul_seismo_depths
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
# get neotectonic domains class and Mmax from zone centroid
###############################################################################
# get path to reference shapefile
shapepath = open('..//reference_shp.txt').read()

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

# loop through Mcomp zones
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
    
    # set completeness based on Leonard 2018 Mcomp model
    if code.startswith('SWA'):
        ycomp.append('1960;1960')
        mcomp.append('3.0;3.0')
        
    elif code.startswith('SA'):
        matchidx = 5
        print 'Fixing index: ', code
        ycomp.append('1966;1966')
        mcomp.append('3.0;3.0')
        
    elif code.startswith('SEA'):
        ycomp.append('1966;1966')
        mcomp.append('2.9;2.9')

    elif code.startswith('WA'):
        ycomp.append('1980;1980')
        mcomp.append('3.1;3.1')
        
    elif code.startswith('EA'):
        ycomp.append('1975;1975')
        mcomp.append('3.0;3.0')
        
    else:
        ycomp.append('1980;1980')
        mcomp.append('3.5;3.5')
            
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
        #bval_fix.append(neo_bval[matchidx])
        #bval_sig_fix.append(bval_sig[matchidx])
        bval_fix.append(-99)
        bval_sig_fix.append(-99)
        #nclass.append(-99)

dep_b = array(dep_b)
    
###############################################################################
# load Rajabi SHMax vectors 
###############################################################################

shmax_pref, shmax_sig = get_aus_shmax_vectors(src_codes, shapes)

# fix preferred upper/lower seismo depths from Domains
usd, lsd = get_ul_seismo_depths(src_codes, usd, lsd)

###############################################################################
# get rate adjustment factors 
###############################################################################

rte_adj_fact = ones_like(usd)
              
###############################################################################
# write initial shapefile
###############################################################################

outshp = 'Mcomp_NSHA18_single.shp'
bval_fix = -99 * ones_like(rte_adj_fact)
bval_sig_fix = -99 * ones_like(rte_adj_fact)

build_source_shape(outshp, shapes, src_names, src_codes, zone_class, \
                   rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                   min_rmag, mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                   shmax_pref, shmax_sig, trt_new, domains, prefCat)

