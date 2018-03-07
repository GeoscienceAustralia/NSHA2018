import shapefile
from shapely.geometry import Polygon, Point
from numpy import ones_like, array
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid
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
src_names = get_field_data(sf, 'NAME', 'str')
src_codes = src_names


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
neo_trt  = get_field_data(dsf, 'TRT', 'str')
neo_dep  = get_field_data(dsf, 'DEP_BEST', 'float')
neo_domains = get_field_data(dsf, 'DOMAIN', 'float')
neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
neo_trt = get_field_data(dsf, 'TRT', 'str')
neo_usd = get_field_data(dsf, 'USD', 'float')
neo_lsd = get_field_data(dsf, 'LSD', 'float')
neo_ow_lsd = get_field_data(dsf, 'OW_LSD', 'float')
neo_stk = get_field_data(dsf, 'PREF_STK', 'float')
neo_dip = get_field_data(dsf, 'PREF_DIP', 'float')
neo_rke = get_field_data(dsf, 'PREF_RKE', 'float')

# get domain polygons
dom_shapes = dsf.shapes()

# set feilds to fill
dom = []
mmax = []
trt = []
nclass = []
dep_b = []
mmax = []
usd = []
lsd = []
ow_lsd = []
stk = []
dip = []
rke = []

# loop through Mcomp zones
for code, poly in zip(src_codes, shapes):
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
    
            
    # set dummy values
    if matchidx == -99:
        dom.append(-99)
        mmax.append(-99)
        trt.append(-99)
        dep_b.append(-99)
        mmax.append(-99)
        usd.append(-99)
        lsd.append(-99)
        ow_lsd.append(-99)
        stk.append(-99)
        dip.append(-99)
        rke.append(-99)
        nclass.append(-99)
    # fill real values
    else:
        dom.append(neo_doms[matchidx])
        mmax.append(neo_mmax[matchidx])
        trt.append(neo_trt[matchidx])
        dep_b.append(neo_dep[matchidx])
        mmax.append(neo_mmax[matchidx])
        usd.append(neo_usd[matchidx])
        lsd.append(neo_lsd[matchidx])
        ow_lsd.append(neo_ow_lsd[matchidx])
        stk.append(neo_stk[matchidx])
        dip.append(neo_dip[matchidx])
        rke.append(neo_rke[matchidx])

lsd = array(lsd)
dom = array(dom)
dom[1] = 2.
###############################################################################
# load 2018 completeness models
###############################################################################

ycomp, mcomp, min_rmag = get_completeness_model(src_codes, shapes, dom, 1)

min_rmag[0] = 3.2

###############################################################################
# load Rajabi SHMax vectors 
###############################################################################

shmax_pref, shmax_sig = get_aus_shmax_vectors(src_codes, shapes)

# fix preferred upper/lower seismo depths from Domains
usd, lsd = get_ul_seismo_depths(src_codes, usd, lsd)

###############################################################################
# make some adjustments
###############################################################################

rte_adj_fact = ones_like(usd)
dep_u = 0. * ones_like(usd)
dep_l = 20. * ones_like(usd)
nclass = dom
              
###############################################################################
# write initial shapefile
###############################################################################

outshp = 'Mcomp_NSHA18_test.shp'
bval_fix = -99 * ones_like(rte_adj_fact)
bval_sig_fix = -99 * ones_like(rte_adj_fact)

build_source_shape(outshp, shapes, src_names, src_codes, nclass, \
                   rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                   min_rmag, mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, stk, dip, rke, \
                   shmax_pref, shmax_sig, trt, dom, prefCat)

