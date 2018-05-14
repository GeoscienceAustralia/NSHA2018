import shapefile
from shapely.geometry import Polygon
from numpy import array, ones_like, nan
try:
    from tools.nsha_tools import get_field_data
    from tools.source_shapefile_builder import get_preferred_catalogue, \
                                               get_completeness_model, get_aus_shmax_vectors, \
                                               get_rate_adjust_factor, build_source_shape, \
                                               aggregate_intraslab_sources
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

###############################################################################

''' START MAIN CODE HERE '''
   
###############################################################################
# parse Domains shp and prep data
###############################################################################

domshp = 'Domains_NSHA18_Merged.shp'

print 'Reading source shapefile...'
sf = shapefile.Reader(domshp)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get field data
src_codes = get_field_data(sf, 'code', 'str')
src_names = get_field_data(sf, 'Name', 'str')
domains = get_field_data(sf, 'DOMAIN', 'float')
mmax = get_field_data(sf, 'max_mag', 'float')
trt = get_field_data(sf, 'trt', 'str')
usd = get_field_data(sf, 'usd', 'float')
lsd = get_field_data(sf, 'lsd', 'float')
hd = get_field_data(sf, 'hd1', 'float')
stk = get_field_data(sf, 'strike1', 'float')
dip = get_field_data(sf, 'dip1', 'float')
rke = get_field_data(sf, 'rake1', 'float')

# set domain for unset domains
trt_new = []
for i in range(0,len(trt)):
    if trt[i] == 'Active':
        domains[i] = 9
    elif trt[i] == 'Extended' and domains[i] == 0.:
        domains[i] = 7
    elif trt[i] == 'Interface':
        domains[i] = 10
    elif trt[i] == 'Intraslab':
        domains[i] = 11
    
    if trt[i] == 'NCratonic':
        trt_new.append('Non_cratonic')
    else:
        trt_new.append(trt[i])

zone_class = list(domains)[:]
# reset Gawler Craton to Flinders due to b-value similarities
zone_class[62] = 2.

# reset W Tasmania to Flinders due to b-value similarities (Delamerian Orogen)
zone_class[61] = 2.

# set pref strike/dip/rake
pref_stk = []
pref_dip = []
pref_rke = []
for i in range(0,len(stk)):
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
        if trt_new[i] == 'Cratonic':
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


###############################################################################
# get preferred catalogues 
###############################################################################

# get preferred catalogues for each zone
prefCat = get_preferred_catalogue(domshp)

# fix catalogue for source zones
prefCat[55] = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
prefCat[56] = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
    
###############################################################################
# load 2018 completeness models
###############################################################################
single_mc = 0
ycomp, mcomp, min_rmag = get_completeness_model(src_codes, shapes, domains, single_mc)
    
# use manual modification
for i in range(0,len(trt_new)):
    if trt_new[i] == 'Active':
        min_rmag[i] = 5.7

min_rmag[3] = 6.1 # TAFS
min_rmag[12] = 6.1 # NBT
min_rmag[11] = 6.0 # NBOT
min_rmag[16] = 5.6 # BNBD
min_rmag[26] = 3.5 # NWO
min_rmag[50] = 3.2 # CARP
min_rmag[51] = 3.5 # EAPM
min_rmag[52] = 3.3 # KMBY
min_rmag[55] = 3.3 # NACR
min_rmag[54] = 3.3 # NAOR
min_rmag[49] = 3.3 # PLBR
min_rmag[53] = 3.5 # WAPM
min_rmag[48] = 3.2 # YLGN
min_rmag[56] = 3.2 # WAEP
min_rmag[66] = 3.5 # NWB1

# SEOB - multi-corner
ycomp[59] = '1980;1964;1900'
mcomp[59] = '3.5;5.0;6.0'
min_rmag[59] = 3.5

###############################################################################
#  set intraslab aggregation class
###############################################################################
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
# load Rajabi SHMax vectors 
###############################################################################

shmax_pref, shmax_sig = get_aus_shmax_vectors(src_codes, shapes)

###############################################################################
# get rate adjustment factors 
###############################################################################

origshp = 'Domains_Sep2011_edit.shp'
newField = 'code'
origField = 'CODE'
rte_adj_fact = get_rate_adjust_factor(domshp, newField, origshp, origField)
              
###############################################################################
# write initial shapefile
###############################################################################
if single_mc == 1:
    outshp = 'Domains_NSHA18_single_Mc.shp'
elif single_mc == 0:
    outshp = 'Domains_NSHA18_multi_Mc.shp'
    
bval_fix = -99 * ones_like(rte_adj_fact)
bval_sig_fix = -99 * ones_like(rte_adj_fact)

build_source_shape(outshp, shapes, src_names, src_codes, zone_class, \
                   rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                   min_rmag, mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                   shmax_pref, shmax_sig, trt_new, domains, prefCat)


