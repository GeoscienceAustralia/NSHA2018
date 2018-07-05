from numpy import array, arange, argsort, where, delete, hstack, sqrt, \
                  unique, mean, percentile, log10, ceil, floor, \
                  nan, isnan, around, diff, interp, exp, ones_like
from os import path, sep, mkdir, getcwd, walk
from shapely.geometry import Point, Polygon
#from osgeo import ogr
from datetime import datetime
from sys import argv
import shapefile
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar, style
from mpl_toolkits.basemap import Basemap
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
from tools.nsha_tools import toYearFraction, get_shapely_centroid
from tools.mfd_tools import * # get_mfds, get_annualised_rates, fit_a_value, parse_hmtk_cat, parse_hmtk_cat
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

style.use('classic')

# import non-standard functions
try:
    from catalogue_tools import weichert_algorithm, aki_maximum_likelihood, bval2beta
    from oq_tools import get_oq_incrementalMFD, beta2bval #, bval2beta
    from mapping_tools import get_field_data, get_field_index, drawoneshapepoly, \
                              drawshapepoly, labelpolygon, get_WGS84_area
    #from catalogue.parsers import parse_ggcat
    from catalogue.writers import ggcat2ascii
    from mag_tools import nsha18_bilin_mw2ml
    from gmt_tools import cpt2colormap 
    from misc_tools import remove_last_cmap_colour, get_log_xy_locs
    
    #from misc_tools import listdir_extension
    #from make_nsha_oq_inputs import write_oq_sourcefile
except:
    cwd = getcwd().split(sep)
    pythonpath = sep.join(pt[0:-3])+sep+'tools'
    print '\nSet environmental variables, e.g.:\n\nexport PYTHONPATH='+pythonpath+':$PYTHONPATH\n'

def timedelta2days_hours_minutes(td):
    return td.days, td.seconds//3600, (td.seconds//60)%60
        
###############################################################################
# set defaults
###############################################################################

bin_width = 0.1
shpfile  = 'shapefiles/Other/Mcomp_NSHA18_multi.shp'

magLabels = ['Original Magnitude', 'Original Magnitude (Revised ML)', \
             'Alternate MW (Bilinear Empirical)', 'Alternate MW (Quadratic Empirical)', \
             'Preferred MW (Quadratic Simulated)']

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_10.cpt'
else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//Paired_10.cpt'
    
ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

###############################################################################
# parse shapefile and make shapely objects
###############################################################################

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
polygonsCopy = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    polygonsCopy.append(Polygon(poly.points))
    
# get input arrays from shapefile
src_code = get_field_data(sf, 'CODE', 'str')
src_name = get_field_data(sf, 'SRC_NAME', 'str')
src_class = get_field_data(sf, 'CLASS', 'str')
src_rte_adj = get_field_data(sf, 'RTE_ADJ_F', 'float')
src_usd = get_field_data(sf, 'USD', 'float')
src_lsd = get_field_data(sf, 'LSD', 'float')
src_overwrite_lsd = get_field_data(sf, 'OW_LSD', 'float')
src_mmin = get_field_data(sf, 'MIN_MAG', 'float')
src_mmin_reg = get_field_data(sf, 'MIN_RMAG', 'float')
src_mmax = get_field_data(sf, 'MMAX_BEST', 'float')
src_mmax_u = get_field_data(sf, 'MMAX_UPPER', 'float')
src_mmax_l = get_field_data(sf, 'MMAX_LOWER', 'float')
src_bval = get_field_data(sf, 'BVAL_BEST', 'float')
src_bval_u = get_field_data(sf, 'BVAL_UPPER', 'float')
src_bval_l = get_field_data(sf, 'BVAL_LOWER', 'float')
src_n0 = get_field_data(sf, 'N0_BEST', 'float')
src_n0_u = get_field_data(sf, 'N0_UPPER', 'float')
src_n0_l = get_field_data(sf, 'N0_LOWER', 'float')
src_bval_fix = get_field_data(sf, 'BVAL_FIX', 'float')
src_bval_fix_sd = get_field_data(sf, 'BVAL_FIX_S', 'float') # too many chars - does not recognise "D"
src_mcomp = get_field_data(sf, 'MCOMP', 'str')
src_ycomp = get_field_data(sf, 'YCOMP', 'str')
src_shmax = get_field_data(sf, 'SHMAX', 'float')
src_shm_sig = get_field_data(sf, 'SHMAX_SIG', 'float')
src_ymax = get_field_data(sf, 'CAT_YMAX', 'float')
src_cat = get_field_data(sf, 'CAT_FILE', 'str')
sortind = argsort(src_code)

###############################################################################
# parse test catalogue
###############################################################################

def parse_alt_mag_catalogue(hmtk_csv):
    from datetime import datetime
    
    lines = open(hmtk_csv).readlines()[1:]
    
    # set arrays
    altMWdict = []
    
    # fill arrays
    for line in lines:
        data = line.strip().split(',')
        dateStr = data[0]
        lon = float(data[9])
        lat = float(data[10])
        dep = float(data[14])
        mx_origML = float(data[19])
        mx_revML = float(data[21])
        mw_pref = float(data[22])
        mw_qds = float(data[25])
        mw_ble = float(data[23])
        mw_qde = float(data[24])
        
        # set empirical alt MW
        if mw_qds == mw_pref:
            mw_alt_ble = mw_ble
            mw_alt_qde = mw_qde
        else:
            mw_alt_ble = mw_pref
            mw_alt_qde = mw_pref
            
        # make datetime object
        try:
            evdt = datetime.strptime(dateStr, '%Y%m%d%H%M')
        except:
            print dateStr
            evdt = datetime.strptime(dateStr, '%Y-%m-%d %H:%M')
        ev_date = evdt
        
        tmpdict = {'datetime':ev_date, 'lon':lon, 'lat':lat, 'dep':dep,
                   'mx_origML':mx_origML, 'mx_revML': mx_revML,
                   'mw_pref':mw_pref, 'mw_qds':mw_qds, 'mw_qde':mw_qde, 
                   'mw_ble':mw_ble, 'mw_alt_ble':mw_alt_ble, 
                   'mw_alt_qde':mw_alt_qde, 'prefmag':mw_pref}
                   	
        altMWdict.append(tmpdict)
        
    return altMWdict, len(altMWdict)

###############################################################################
# initiate new arrays for writing new shpfile
###############################################################################

new_bval_b = src_bval  
new_bval_l = src_bval_l
new_bval_u = src_bval_u
new_n0_b = src_n0
new_n0_l = src_n0_l
new_n0_u = src_n0_u

# reset Mmin to 4.8
#print '!!!Setting Mmin = 4.5!!!'
#src_mmin = 4.5 * ones_like(src_mmin)
#src_mmin_reg = 4. * ones_like(src_mmin_reg)

# set arrays for testing
bval_vect = []
bsig_vect = []

srcidx = range(len(src_code))

###############################################################################
# parse NSHA-Cat and ISC-GEM catalogues
###############################################################################

# parse NSHA-Cat catalogue
hmtk_csv = path.join('..','..','catalogue','data','NSHA18CAT_V0.2_hmtk_declustered.csv')
nshaCat, neq = parse_alt_mag_catalogue(hmtk_csv)
nshaMaxYear = toYearFraction(nshaCat[-1]['datetime'])

###############################################################################
# loop through source zones
###############################################################################
mfdtxt = 'SRCZONE,SRCAREA,MAGTYPE,NEQ,A0,BVAL,BVALSIG\n'

src_area = [] 
fig = plt.figure(1, figsize=(16, 10))
k = 0
for i in srcidx:
    print '\nFitting MFD for', src_code[i]
    
    if i == 2 or i == 3:
        k += 1
        ax = fig.add_subplot(1,2,k)
    
    # get completeness periods for zone
    ycomps = array([int(x) for x in src_ycomp[i].split(';')])
    mcomps_mw = array([float(x) for x in src_mcomp[i].split(';')])
    mcomps_mw[0] = 2.95
    mcomps_ml = [nsha18_bilin_mw2ml(mw) for mw in mcomps_mw]
    mcomps_ml[0] = 2.95
    
    print 'if ml based, adjust Mcomp'
    
    # get mag range for zonea
    mcompmin_mw = min(mcomps_mw)
    mcompmin_ml = min(mcomps_ml)
            
    # convert min Mx to MW
    mcompminml = around(ceil(mcompmin_ml*10.) / 10., decimals=1)
    mcompminml = 2.95
    mcompminmw = around(ceil(mcompmin_mw*10.) / 10., decimals=1)
    mcompminmw = 2.95
    mrng_mw = arange(mcompminmw-bin_width/2, src_mmax[i], bin_width)
    mrng_ml = arange(mcompminml-bin_width/2, src_mmax[i], bin_width)
            
    # set null values to avoid plotting issues later
    bval = 1.
    bval_sig = 0.1
    new_bval_b[i] = 1.0
    new_n0_b[i] = 1E-30
    
    # set beta params       
    beta = bval2beta(bval)
    sigbeta = bval2beta(bval_sig)
    
    # set polygons  
    poly = polygons[i]
    
    # get area (in km**2) of sources for normalisation
    src_area.append(get_WGS84_area(poly))
    
    ###############################################################################
    # set preferred catalogue for each source
    ###############################################################################
    
    # mw alt based om mw_ble; pref_mw based on mw_qds
    magKeys = ['mx_origML', 'mx_revML', 'mw_alt_ble', 'mw_alt_qde', 'mw_pref']
    
    # loop through mag keys and get events
    for j, mk in enumerate(magKeys):
        for n in range(0, len(nshaCat)):
            nshaCat[n]['prefmag'] = nshaCat[n][mk]
            
        # set Mcomps
        if mk == 'mx_origML' or mk == 'mx_revML':
            mcomps = mcomps_ml
            mrng = mrng_ml
        else:
            mcomps = mcomps_mw
            mrng = mrng_mw
                
        # now get events within zone of interest
        mvect, mxvect, tvect, dec_tvect, ev_dict = get_events_in_poly_simple(nshaCat, poly)
        print 'NEQ Before =', len(mvect)
            
        # remove incomplete events based on new MW estimates (mvect)
        mvect, mxvect, tvect, dec_tvect, ev_dict, out_idx, ev_out = \
             remove_incomplete_events(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, bin_width)
            
        # check to see if mvect still non-zero length after removing incomplete events
        print 'NEQ After =', len(mvect)
        print '\n'+mk
        
        # get bval for combined zones data - uses new MW estimates ("total_mvect") to do cleaning
        bval, beta, sigb, sigbeta, fn0, cum_rates, ev_out, err_up, err_lo = \
              get_mfds(mvect, mxvect, tvect, dec_tvect, ev_dict, \
                       mcomps, ycomps, nshaMaxYear, mrng, src_mmax[i], src_mmin_reg[i], \
                       src_bval_fix[i], src_bval_fix_sd[i], bin_width, poly)
        
        # mk txt:'SRCZONE,SRCAREA,MAGTYPE,NEQ,A0,BVAL,BVALSIG\n'
        mfdtxt += ','.join((src_code[i], str(round(src_area[i])), mk, str(len(mvect)), \
                            str('%0.2f' % log10(fn0)), str('%0.2f' % bval), str('%0.3f' % sigb))) + '\n'
        ###############################################################################
        # start making plots
        ###############################################################################
        
        if i == 2 or i == 3:
            # now plt unique values for current source
            uidx = unique(cum_rates[::-1], return_index=True, return_inverse=True)[1]
                    
            plt.semilogy(mrng[::-1][uidx], cum_rates[::-1][uidx], 'o', c=cs[2*j], ms=7, mec=cs[2*j])
                    
            # get betacurve for source
            mpltmin = round(mcomps[0],1) - bin_width/2.
            betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mpltmin, mrng[-1], bin_width)
            
            plt.semilogy(mfd_mrng[:-1], betacurve[:-1], '-', lw=1.5, c=cs[2*j+1], label=magLabels[j], zorder=1000)
    
    ###############################################################################
    # finish mfd
    ###############################################################################
    if i == 2 or i == 3:
        if i == 2:
            plt.ylabel('Cumulative Rate (/yr)', fontsize=16)
            plt.text(2.1, 85, 'a)', fontsize=15, va='top', ha='left')
        else:
            plt.text(2.1, 85, 'b)', fontsize=15, va='top', ha='left')
        plt.xlabel('Magnitude (MW)', fontsize=16)
        
        plt.ylim([2.5, 7])
        plt.ylim([1E-3, 1E2])
        plt.grid(which='both')
        plt.legend(loc=1, fontsize=13)
    
        ###############################################################################
        # make map inset
        ###############################################################################
        
        #a = plt.axes([.05, .05, .25, .25])
        
        axins = inset_axes(ax,
                       width="39%",  # width = 30% of parent_bbox
                       height=2.1,  # height : 1 inch
                       loc=3)
        
        m = Basemap(projection='merc',\
                    llcrnrlon=111,llcrnrlat=-45, \
                    urcrnrlon=156,urcrnrlat=-9,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
                    
        m.drawmapboundary(fill_color='0.8', zorder=0)
        m.fillcontinents(color='w', lake_color='0.8', zorder=1)
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        
        # fill main area
        drawoneshapepoly(m, plt, sf, 'CODE', src_code[i], lw=1.5, col='r')
        
        # save figure
#plt.savefig(path.join('cat_mfd_test', src_code[i]+'_mfdplt.png'), fmt='png', bbox_inches='tight')
plt.savefig(path.join('cat_mfd_test', 'mfdplt.png'), fmt='png', bbox_inches='tight')
    
###############################################################################
# write csv
############################################################################### 
csvfile = path.join('cat_mfd_test', 'cat_mfd_test.csv')
f = open(csvfile, 'wb')
f.write(mfdtxt)
f.close()

plt.show()                
