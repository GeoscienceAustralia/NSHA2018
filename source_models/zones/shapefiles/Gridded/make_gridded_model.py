import shapefile
from shapely.geometry import Polygon, Point
from numpy import array, ones_like, nan, linspace, arange
from sys import argv

from tools.nsha_tools import get_field_data
from tools.source_shapefile_builder import get_preferred_catalogue, \
                                           get_completeness_model, get_aus_shmax_vectors, \
                                           get_rate_adjust_factor, build_source_shape, \
                                           get_ul_seismo_depths, get_simple_neotectonic_domain_params, \
                                           aggregate_intraslab_sources

'''
except:
    print('Add PYTHONPATH to NSHA18 root directory')
'''
###############################################################################

''' START MAIN CODE HERE '''
   
###############################################################################
# make gridded polygons
###############################################################################
res = int(argv[1]) # degrees
r2 = res / 2.
outshp = 'gridded_polygons_'+str(res)+'deg.shp'

print('Making grids shapefile...')

bbox = '110.0/156.0/-45.0/-9.0' # map boundary - lon1/lon2/lat1/lat2
bbox = bbox.split('/')
minlon = float(bbox[0])
maxlon = float(bbox[1])
minlat = float(bbox[2])
maxlat = float(bbox[3])

# make first grid
xrng = arange(minlon-r2, maxlon+r2, res)
yrng = arange(minlat-r2, maxlat+r2, res) # check 1st + shouldn't be a -

polygons = []
zcode = []
r2 = res / 2.
for x in xrng:
    for y in yrng:
        #print('\n'+str(x)+' '+str(y))
        
        # make points
        p1 = [x-r2, y-r2]
        p2 = [x+r2, y-r2]
        p3 = [x+r2, y+r2]
        p4 = [x-r2, y+r2]
        
        pointList = [p1, p2, p3, p4, p1]
        #polygons.append(Polygon(pointList))
        polygons.append(pointList)
        
        zcode.append(str(x)+'E_'+str(abs(y))+'S')

# make offset grid
xrng = arange(minlon, maxlon+r2, res)
yrng = arange(minlat, maxlat+r2, res)

for x in xrng:
    for y in yrng:
        #print('\n'+str(x)+' '+str(y))
        
        # make points
        p1 = [x-r2, y-r2]
        p2 = [x+r2, y-r2]
        p3 = [x+r2, y+r2]
        p4 = [x-r2, y+r2]
        
        pointList = [p1, p2, p3, p4, p1]
        #polygons.append(Polygon(pointList))
        polygons.append(pointList)
        
        zcode.append(str(x)+'E_'+str(abs(y))+'S')
        
#w = shapefile.Writer(shapefile.POLYGON)
w = shapefile.Writer(outshp[:-4], shapeType=5)
w.field('SRC_NAME','C','12')
w.field('CODE','C','12')

for i, poly in enumerate(polygons):
    
        # set shape polygon
        #w.line(parts=[poly], shapeType=shapefile.POLYGON)
        w.poly([poly])
            
        # write new records
        if i >= 0:
            w.record(zcode[i], zcode[i])

# now save area shapefile
#w.save(outshp)
w.close()

###############################################################################
# read new shapefile
###############################################################################

print('Reading source shapefile...')
sf = shapefile.Reader(outshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get field data
src_codes = get_field_data(sf, 'CODE', 'str')
src_names = get_field_data(sf, 'SRC_NAME', 'str')      

###############################################################################
# load neotectonic domains parameters
###############################################################################

# set domestic domain numbers based on neotectonic domains
refshpfile = '..//reference_shp.txt'
neo_domains, neo_min_rmag, neo_mmax, neo_trt, neo_bval_fix, neo_bval_sig_fix, neo_usd, neo_lsd, neo_dep_b, neo_dep_u, neo_dep_l \
    = get_simple_neotectonic_domain_params(sf, refshpfile)
zone_class = list(neo_domains)[:]

'''
# reset Gawler Craton to Flinders due to b-value similarities
zone_class[4] = 2.

# reset West Coast to extended
zone_class[14] = 6.
domains[14] = 6.
mmax[14] = 7.7

# reset West Coast Passive Margin to extended
zone_class[15] = 7.
domains[15] = 7.

# reset Ottway/Gippsland to extended
zone_class[0] = 5.
domains[0] = 5.

# reset Southern Oceanic buffer
zone_class[19] = 8.
domains[19] = 8.

# reset Tasmania buffer
zone_class[11] = 2.
domains[11] = 2.


###############################################################################
#  set intraslab aggregation class
###############################################################################
print('\n!!! REMEMBER TO RESET SRM_200_300 SOURCE CODE !!!!\n'
new_src_codes = []
for i, src_code in enumerate(src_codes):
    
    # fix source codes on the fly
    if src_code.startswith('TMR') or src_code.startswith('SRM'):
        if src_code[3] != '_':
            src_code = src_code[:3]+'_'+src_code[3:]
    
    new_src_codes.append(src_code)
    
    # now match zone class based on lookup table
    zone_class[i] =  aggregate_intraslab_sources(src_code, zone_class[i])

# supplant new source codes
del src_codes
src_codes = array(new_src_codes).copy()

###############################################################################
#  set pref strike/dip/rake
###############################################################################

pref_stk = []
pref_dip = []
pref_rke = []
    if stk[i] == 0.0 and domains[i] <= 8:
        pref_stk.append(-999)
        pref_dip.append(-999)
        pref_rke.append(-999)
    else:
        pref_stk.append(stk[i])
        pref_dip.append(dip[i])
        pref_rke.append(rke[i])
        
# set depth parameters
dep_b = []
dep_u = []
dep_l = []
for i in range(0,len(lsd)):
    if domains[i] <= 8:
        lsd[i] = 20.
        if trt[i] == 'Cratonic':
            dep_b.append(5.0)
            dep_u.append(2.5)
            dep_l.append(10.)
        else:
            dep_b.append(10.)
            dep_u.append(5.0)
            dep_l.append(15.)
    
    # don't care about depth logic tree outside Australia
    else:
        dep_b.append(hd[i])
        dep_u.append(-999)
        dep_l.append(-999)

# fix preferred upper/lower seismo depths from Domains
usd, lsd = get_ul_seismo_depths(src_codes, usd, lsd)
'''

###############################################################################
# get preferred catalogues 
###############################################################################

# get preferred catalogues for each zone
prefCat = get_preferred_catalogue(outshp)
'''
# fix catalogue for source zones
prefCat[44] = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
prefCat[2] = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
'''
###############################################################################
# load 2018 completeness models
###############################################################################
single_mc = 0
ycomp, mcomp, min_rmag_ignore = get_completeness_model(src_codes, shapes, neo_domains, single_mc)
#ycomp, mcomp, min_rmag_ignore = get_completeness_model_point(src_codes, shapes, single_mc)

# use values from Domains model instead
min_rmag = neo_min_rmag

# use manual modification
for i in range(0,len(neo_trt)):
    if neo_trt[i] == 'Active':
        min_rmag[i] = 5.75
    elif neo_trt[i] == 'Intraslab':
        min_rmag[i] = 5.75
    elif neo_trt[i] == '':
        min_rmag[i] = 5.75
'''
min_rmag[23] = 6.1 # TAFS
min_rmag[31] = 6.0 # NBOT
min_rmag[32] = 6.1 # NBT
min_rmag[44] = 3.8 # TP
min_rmag[15] = 3.5 # ZN7d
min_rmag[0] = 3.0 # ZN7d
min_rmag[14] = 3.2 # ZN6b


#min_rmag[46] = 3.8 # NWO
min_rmag[45] = 3.5 # NECS
#min_rmag[50] = 3.2 # CARP
min_rmag[19] = 3.5 # SEOB
min_rmag[18] = 3.5 # SWOB



# SEOB - multi-corner
ycomp[19] = '1980;1964;1900'
mcomp[19] = '3.5;5.0;6.0'

ycomp[18] = '1980;1964;1900'
mcomp[18] = '3.5;5.0;6.0'

ycomp[46] = '1980;1964;1900'
mcomp[46] = '3.5;5.0;6.0'

ycomp[13] = '1980;1964;1900'
mcomp[13] = '3.5;5.0;6.0'
'''


###############################################################################
# load Rajabi SHMax vectors 
###############################################################################

shmax_pref, shmax_sig = get_aus_shmax_vectors(src_codes, shapes)

###############################################################################
# get rate adjustment factors 
###############################################################################
'''
origshp = 'ARUP_NSHA18_FIXEDSHAPES.shp'
newField = 'CODE'
origField = 'CODE'
rte_adj_fact = get_rate_adjust_factor(domshp, newField, origshp, origField)
'''           
###############################################################################
# write initial shapefile
###############################################################################

#outshp = 'ARUP_NSHA18.shp'
bval_fix = -99 * ones_like(array(min_rmag))
bval_sig_fix = -99 * ones_like(array(min_rmag))
rte_adj_fact = ones_like(array(min_rmag))
pref_stk = -99 * ones_like(array(min_rmag))
pref_dip = -99 * ones_like(array(min_rmag))
pref_rke = -99 * ones_like(array(min_rmag))

build_source_shape(outshp, shapes, src_names, src_codes, zone_class, \
                   rte_adj_fact, neo_dep_b, neo_dep_u, neo_dep_l, neo_usd, array(neo_lsd), \
                   min_rmag, neo_mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                   shmax_pref, shmax_sig, neo_trt, neo_domains, prefCat)


# write projection file
print(outshp)
prjfile = outshp.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()
