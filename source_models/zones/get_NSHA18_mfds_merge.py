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
from matplotlib import colors, colorbar
from mpl_toolkits.basemap import Basemap
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
from tools.nsha_tools import toYearFraction, get_shapely_centroid
from mfd_tools import * # get_mfds, get_annualised_rates, fit_a_value, parse_hmtk_cat

# import non-standard functions
try:
    from catalogue_tools import weichert_algorithm, aki_maximum_likelihood, bval2beta
    from oq_tools import get_oq_incrementalMFD, beta2bval#, bval2beta
    from mapping_tools import get_field_data, get_field_index, drawoneshapepoly, \
                              drawshapepoly, labelpolygon, get_WGS84_area
    #from catalogue.parsers import parse_ggcat
    from catalogue.writers import ggcat2ascii
    
    
    #from misc_tools import listdir_extension
    #from make_nsha_oq_inputs import write_oq_sourcefile
except:
    cwd = getcwd().split(sep)
    pythonpath = sep.join(pt[0:-3])+sep+'tools'
    print '\nSet environmental variables, e.g.:\n\nexport PYTHONPATH='+pythonpath+':$PYTHONPATH\n'

def timedelta2days_hours_minutes(td):
    return td.days, td.seconds//3600, (td.seconds//60)%60
        
###############################################################################
# do check to see if running single source only
###############################################################################

#try:
# parse param file
paramfile = argv[1]

try:
    single_zone = array([argv[2].upper()])
    single_src = True
except:
    single_src = False

#print '!!!Temporary test - setting weights of upper and lower curves to zero!!!'
#bestcurve = True
###############################################################################
# get parameter values
###############################################################################

# load param file
lines = open(paramfile).readlines()
rootfolder  = lines[0].split('=')[-1].strip()
hmtk_csv    = lines[1].split('=')[-1].strip()
dec_flag    = lines[2].split('=')[-1].strip() # decluster flag
shpfile     = lines[3].split('=')[-1].strip()
outfolder   = path.join(rootfolder, lines[4].split('=')[-1].strip())
outsrcshp   = lines[5].split('=')[-1].strip()
bin_width   = float(lines[6].split('=')[-1].strip())
#bval_shp    = lines[7].split('=')[-1].strip()

# get export folder
now = datetime.now()
shpname = path.split(shpfile)
'''
if dec_flag == 'True':
    #outfolder = '_'.join((shpname[-1][0:-4], 'DEC', now.strftime('%Y-%m-%d')))
    outfolder += '_'.join(('DEC', now.strftime('%Y-%m-%d')))
else:
    outfolder += '_' + now.strftime('%Y-%m-%d')
'''
# check to see if exists
if path.isdir(outfolder) == False:
    mkdir(outfolder)

shpfolder = path.split(shpfile)

###############################################################################
# parse shapefile and make shapely objects
###############################################################################

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get input arrays from shapefile
src_code = get_field_data(sf, 'CODE', 'str')
src_name = get_field_data(sf, 'SRC_NAME', 'str')
src_class = get_field_data(sf, 'CLASS', 'str')
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
src_ymax = get_field_data(sf, 'YMAX', 'float')
src_cat = get_field_data(sf, 'CAT_FILE', 'str')
sortind = argsort(src_code)

# initiate new arrays for writing new shpfile
new_bval_b = src_bval  
new_bval_l = src_bval_l
new_bval_u = src_bval_u
new_n0_b = src_n0
new_n0_l = src_n0_l
new_n0_u = src_n0_u

# reset Mmin to 4.8
print '!!!Setting Mmin = 4.5!!!'
src_mmin = 4.5 * ones_like(src_mmin)
#src_mmin_reg = 4. * ones_like(src_mmin_reg)

# set all plausible depths
depmin = -99.
depmax = 99. # only get events GE 10 km

# set shallow depth range
#depmin = -99.
#depmax = 99. # only get events GE 10 km

# set arrays for testing
bval_vect = []
bsig_vect = []

# if single source remove unnecessary data
if single_src == True:
    srcidx = where(array(src_code) == single_zone)[0]
    ssi = where(sortind == srcidx)[0]

else:
    srcidx = range(len(src_code))

###############################################################################
# parse catalogue
###############################################################################
'''Used to parse GGCat csv - now parse HMTK csv'''

ggcat, neq = parse_hmtk_cat(hmtk_csv)
    
# get max decimal year and round up!
lastRec = ggcat[-1]
year_max = lastRec['year'] + lastRec['month']/12.

# apply Hadi's ML-MW conversion to mcomps
def convert_mcomps(mcomps):
    '''
    a1 = 0.66199378
    a2 = 1.2156352
    a3 = 1.2156352
    mx = 4.5
    my = a1 * mx + a2
    
    idx = where(mcomps <= mx)[0]
    mcomps[idx] = a1 * mcomps[idx] + a2
    idx = where(mcomps > mx)[0]
    mcomps[idx] = a3 * (mcomps[idx] - mx) + my
    '''
    return mcomps

###############################################################################
# get unique zone classes and loop through to merge zones of similar class 
###############################################################################

unique_classes = unique(array(src_class))
class_bval = []
class_bval_sig = []
class_cum_rates = []
class_mrng = []
class_err_up = []
class_err_lo = []
class_idxs = []
class_codes = []
class_mvect = []
class_lavect = []
class_lovect = []
class_fn0 = []
class_area = []

for uclass in unique_classes:
    
    total_mvect = []
    total_mxvect = []
    total_tvect = []
    total_dec_tvect = []
    total_ev_dict = []
    #out_idx = []
    class_idx = []
    class_code = []
    lavect = []
    lovect = []
    L08_b = False
    Aki_ML = False
    Weichert = False
    mcomps = [-9999]
    cum_area = 0
    class_mmin_reg = 9.0
    class_mmax = nan
    
    print '\nCalculating b-value for class:', uclass
            
    
    ###############################################################################
    # loop thru zones 
    ###############################################################################
    
    # loop thru source zones
    for i in srcidx:
        
        if src_class[i] == uclass:
            class_idx.append(i)
            class_code.append(src_code[i])
        
            # get polygon of interest
            poly = polygons[i]
            
            # get cumulative class area
            cum_area += get_WGS84_area(poly)
            
            print '    Compiling data from:', src_code[i]
            
            # definitions of arrays extracted from catalogue
            '''
            mvect: preferred MW
            mxvect: preferred original magnitudes
            tvect: datetime array
            dec_tvect: decimal datetime
            ev_dict: event dictionary
            '''
            
            # get earthquakes within source zones
            mvect, mxvect, tvect, dec_tvect, ev_dict \
                   = get_events_in_poly(ggcat, poly, depmin, depmax)
            
            # stack records into total arrays
            total_mvect = hstack((total_mvect, mvect))
            total_mxvect = hstack((total_mxvect, mxvect))
            total_tvect = hstack((total_tvect, tvect))
            total_dec_tvect = hstack((total_dec_tvect, dec_tvect))
            total_ev_dict = hstack((total_ev_dict, ev_dict))
                    
            ###############################################################################
            # get earthquakes that pass completeness in merged zones
            ###############################################################################
            
            # get most conservative completeness for given geologic class
            temp_mcomp = array([float(x) for x in src_mcomp[i].split(';')])
            if temp_mcomp[-1] > mcomps[0]:
                ycomps = array([int(x) for x in src_ycomp[i].split(';')])
                mcomps = array([float(x) for x in src_mcomp[i].split(';')])
            
            # get mag range for mfd
            mcompmin = min(mcomps)
            
            # convert min Mx to MW
            #mcompminmw = aus_ml2mw(mcompmin) - not needed now that mcomps in MW
            
            mcompminmw = around(ceil(mcompmin*10.) / 10., decimals=1)
            mrng = arange(mcompminmw-bin_width/2, src_mmax[i], bin_width)
            
            # remove events with NaN mags
            didx = where(isnan(total_mvect))[0]
            total_tvect = delete(total_tvect, didx)
            total_mvect = delete(total_mvect, didx)
            total_mxvect = delete(total_mxvect, didx)
            total_dec_tvect = delete(total_dec_tvect, didx)
            total_ev_dict = delete(total_ev_dict, didx)
            
            # remove events with M < min mcomps
            didx = where(total_mvect < min(mcomps)-bin_width/2)[0]
            total_tvect = delete(total_tvect, didx)
            total_mvect = delete(total_mvect, didx)
            total_mxvect = delete(total_mxvect, didx)
            total_dec_tvect = delete(total_dec_tvect, didx)
            total_ev_dict = delete(total_ev_dict, didx)
            
            # set min regression magnitude
            if src_mmin_reg[i] < class_mmin_reg:
                class_mmin_reg = src_mmin_reg[i]
                
            # set class mmax
            class_mmax = src_mmax[i]
            
            # set class b-values
            fixed_bval = src_bval_fix[i]
            fixed_bval_sig = src_bval_fix_sd[i]
            
       
    ###############################################################################
    # get b-values from joined zones
    ###############################################################################
    
    # keep original vectors for plotting
    class_orig_mvect = total_mvect
    class_orig_mxvect = total_mxvect
    class_orig_tvect = total_tvect
    class_orig_dec_tvect = dec_tvect
    
    # get bval for combined zones data - uses new MW estimates ("total_mvect") to do cleaning
    bval, beta, sigb, sigbeta, fn0, cum_rates, ev_out, err_up, err_lo = \
          get_mfds(total_mvect, total_mxvect, total_tvect, total_dec_tvect, total_ev_dict, \
                   mcomps, ycomps, year_max, mrng, class_mmax, class_mmin_reg, \
                   fixed_bval, fixed_bval_sig, bin_width, poly)
    
    # get a-value using fixed region class b-value if assigned - need to do this to fit the class rates!
    if not fixed_bval == -99.0:
        
        # remove incomplete events based on original preferred magnitudes (mxvect)
        total_mvect, total_mxvect, total_tvect, total_dec_tvect, total_ev_dict, out_idx, ev_out = \
             remove_incomplete_events(total_mvect, total_mxvect, total_tvect, total_dec_tvect, total_ev_dict, mcomps, ycomps, bin_width)        
        
        # get annualised rates based on preferred MW (mvect)
        cum_rates, cum_num, bin_rates, n_obs, n_yrs = \
            get_annualised_rates(mcomps, ycomps, total_mvect, mrng, bin_width, year_max)
            
        # get index of min reg mag and valid mag bins
        diff_cum = abs(hstack((diff(cum_rates), 0.)))
        midx = where((mrng >= class_mmin_reg-bin_width/2) & (diff_cum > 0.))[0]

        # check if length of midx = 0 and get highest non-zero mag
        if len(midx) == 0:
            midx = [where(isfinite(diff_cum))[0][-1]]
        
        # make sure there is at least 4 observations for a-value calculations
        if len(midx) < 5:
            idxstart = midx[0] - 1
            
            while idxstart >= 0 and len(midx) < 5:
                # if num observations greater than zero, add to midx
                if n_obs[idxstart] > 0:
                    midx = hstack((idxstart, midx))
                    print '    get lower mag M', midx
                    
                idxstart -= 1
        
        # reset fn0 based on fixed b-value        
        fn0 = fit_a_value(fixed_bval, mrng, cum_rates, class_mmax, bin_width, midx)
        
    # add to class arrays - used for plotting later
    class_bval.append(bval)
    class_bval_sig.append(sigb)
    class_cum_rates.append(cum_rates)
    class_mrng.append(mrng)
    class_err_up.append(err_up)
    class_err_lo.append(err_lo)
    class_idxs.append(class_idx)
    class_codes.append(class_code)
    class_mvect.append(total_mvect)
    class_fn0.append(fn0)
    class_area.append(cum_area)
    
    # get event coords for plotting
    for e in total_ev_dict:
        lavect.append(e['lat'])
        lovect.append(e['lon'])
    
    class_lavect.append(lavect)
    class_lovect.append(lovect)
    
###############################################################################
# loops through individual sources using class b-value
###############################################################################

src_area = [] 

for i in srcidx:
    print '\nFitting MFD for', src_code[i]
    
    # get completeness periods for zone
    ycomps = array([int(x) for x in src_ycomp[i].split(';')])
    mcomps = array([float(x) for x in src_mcomp[i].split(';')])
    
    # get mag range for zonea
    mcompmin = min(mcomps)
            
    # convert min Mx to MW
    #mcompminmw = aus_ml2mw(mcompmin) # now defining Mc in terms of MW
    mcompminmw = around(ceil(mcompmin*10.) / 10., decimals=1)
    mrng = arange(mcompminmw-bin_width/2, src_mmax[i], bin_width)
            
    for uc in range(0, len(unique_classes)):
        
        # set b-value and class data
        if src_class[i] == unique_classes[uc]:
            bval = class_bval[uc]
            bval_sig = class_bval_sig[uc]
            
            source_class = unique_classes[uc]
            class_idx = uc
            
            bval_vect.append(bval)
            bsig_vect.append(bval_sig)
            
    # set null values to avoid plotting issues later
    try:
        new_bval_b[i] = bval
    except:
        new_bval_b[i] = 1.0
    
    new_n0_b[i] = 1E-30
    
    # set beta params       
    beta = bval2beta(bval)
    sigbeta = bval2beta(bval_sig)
            
    # get polygon of interest
    poly = polygons[i]
    
    # get area (in km**2) of sources for normalisation
    src_area.append(get_WGS84_area(poly))
    
    # now get events within zone of interest
    mvect, mxvect, tvect, dec_tvect, ev_dict \
        = get_events_in_poly(ggcat, polygons[i], depmin, depmax)
        
    # skip zone if no events pass completeness
    if len(mvect) != 0:
        # assume Banda Sea sources
        if src_mmax[i] == -99:
            src_mmax[i] = 8.0
            src_mmax_l[i] = 7.8
            src_mmax_u[i] = 8.2
    
        # preserve original arrays for plotting
        orig_mvect = mvect
        orig_mxvect = mxvect
        orig_tvect = tvect
        orig_dec_tvect = dec_tvect
        
        # remove incomplete events based on new MW estimates (mvect)
        mvect, mxvect, tvect, dec_tvect, ev_dict, out_idx, ev_out = \
             remove_incomplete_events(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, bin_width)
        
    # check to see if mvect still non-zero length after removing incomplete events
    if len(mvect) != 0:
        
        # get annualised rates based on preferred MW (mvect)
        cum_rates, cum_num, bin_rates, n_obs, n_yrs = \
            get_annualised_rates(mcomps, ycomps, mvect, mrng, bin_width, year_max)
            
        # get index of min reg mag and valid mag bins
        diff_cum = abs(hstack((diff(cum_rates), 0.)))
        midx = where((mrng >= src_mmin_reg[i]-bin_width/2) & (diff_cum > 0.))[0]
        
        # check if length of midx = 0 and get highest non-zero mag
        if len(midx) == 0:
            midx = [where(isfinite(diff_cum))[0][-1]]
        
        # make sure there is at least 4 observations for a-value calculations
        if len(midx) < 5:
            idxstart = midx[0] - 1
            
            while idxstart >= 0 and len(midx) < 5:
                # if num observations greater than zero, add to midx
                if n_obs[idxstart] > 0:
                    midx = hstack((idxstart, midx))
                    print '    get lower mag M2', midx
                    
                idxstart -= 1
        
        # get a-value using region class b-value
        fn0 = fit_a_value(bval, mrng, cum_rates, src_mmax[i], bin_width, midx)
        
        # get zone confidence limits
        err_up, err_lo = get_confidence_intervals(n_obs, cum_rates)
        
        ###############################################################################
        # get upper and lower MFD bounds
        ###############################################################################
        sigbeta173 = 1.73 * sigbeta
        sigb173 = 1.73 * bval_sig
    
        # preallocate data
        N0_lo173 = nan
        N0_up173 = nan
        N0_lo100 = nan
        N0_up100 = nan
    
        if not isnan(bval):
            
            mpltmin = mrng[midx][0]
            dummyN0 = 1.
            
            # get lower + 1 std
            bc_tmp, bc_mrng_lo = get_oq_incrementalMFD(beta+sigbeta, dummyN0, mpltmin, src_mmax_l[i], bin_width)
            # fit to err_lo
            bc_lo100 = (cum_rates[midx][0] - err_lo[midx][0]) * (bc_tmp / bc_tmp[0])
            # solve for N0
            N0_lo100 = 10**(log10(bc_lo100[0]) + beta2bval(beta+sigbeta)*bc_mrng_lo[0])
            
            # get lower + 1.73 std
            bc_tmp, bc_mrng_lo = get_oq_incrementalMFD(beta+sigbeta173, dummyN0, mpltmin, src_mmax_l[i], bin_width)
            # fit to err_lo
            bc_lo173 = (cum_rates[midx][0] - err_lo[midx][0]) * (bc_tmp / bc_tmp[0])
            # solve for N0
            N0_lo173 = 10**(log10(bc_lo173[0]) + beta2bval(beta+sigbeta173)*bc_mrng_lo[0])
            
            # get upper - 1 std
            bc_tmp, bc_mrng_up = get_oq_incrementalMFD(beta-sigbeta, dummyN0, mpltmin, src_mmax_u[i], bin_width)
            # fit to err_up
            bc_up100 = (cum_rates[midx][0] + err_up[midx][0]) * (bc_tmp / bc_tmp[0])
            # solve for N0
            N0_up100 = 10**(log10(bc_up100[0]) + beta2bval(beta-sigbeta)*bc_mrng_up[0])
            
            # get upper - 1.73 std
            bc_tmp, bc_mrng_up = get_oq_incrementalMFD(beta-sigbeta173, dummyN0, mpltmin, src_mmax_u[i], bin_width)
            # fit to err_up
            bc_up173 = (cum_rates[midx][0] + err_up[midx][0]) * (bc_tmp / bc_tmp[0])
            # solve for N0
            N0_up173 = 10**(log10(bc_up173[0]) + beta2bval(beta-sigbeta173)*bc_mrng_up[0])
        
        ###############################################################################
        # fill new values
        ###############################################################################
    
        print 'Filling new values for', src_code[i]
        new_bval_b[i] = bval
        #new_bval_l[i] = bval-sigb173
        #new_bval_u[i] = bval+sigb173
        
        # Use +/- 1 sigma
        new_bval_l[i] = bval+bval_sig # lower curve, so higher b
        new_bval_u[i] = bval-bval_sig # upper curve, so lower b
        
        new_n0_b[i]   = fn0
        new_n0_l[i]   = N0_lo100
        new_n0_u[i]   = N0_up100
        
        ###############################################################################
        # plot earthquakes that pass completeness
        ###############################################################################
        plt.clf()
        fig = plt.figure(i+2, figsize=(20, 9)) # not sure why, but 2nd plot always a dud, so do i+2
        
        # plot original data
        ax = plt.subplot(231)
    
        # cannot plot datetime prior to 1900 - not sure why!
        # plot all events
        dcut = datetime(1900,1,1,0,0)
        didx = where(orig_tvect > dcut)[0]
        h1 = plt.plot(orig_dec_tvect[didx], orig_mvect[didx], 'bo')
        
        # replot those that pass completeness
        didx = where(tvect > dcut)[0]
        h2 = plt.plot(dec_tvect[didx], mvect[didx], 'ro')
        plt.ylabel('Moment Magnitude (MW)')
        plt.xlabel('Date')
        #plt.title(src_code[i] + ' Catalogue Completeness')
        
        
        # now loop thru completeness years and mags
        for yi in range(0, len(ycomps)-1):
            # convert y to datetime
            ydt = datetime(ycomps[yi], 1, 1)
            
            # plt H completeness ranges
            if yi == 0:
                plt.plot([toYearFraction(ydt), 2020], \
                         [mcomps[yi], mcomps[yi]], 'g-', lw=1.5)
            else:
                ydtp = datetime(ycomps[yi-1], 1, 1)
                plt.plot([toYearFraction(ydtp), toYearFraction(ydt)], \
                         [mcomps[yi], mcomps[yi]], 'g-', lw=1.5)
            
            # plt V completeness ranges
            plt.plot([toYearFraction(ydt), toYearFraction(ydt)], \
                     [mcomps[yi], mcomps[yi+1]], 'g-', lw=1.5)
        
        # plt last H completeness range
        ydt = datetime(ycomps[-1], 1, 1)
        if ydt < dcut:
            ydt = dcut
        ydtp = datetime(ycomps[-2], 1, 1)
        plt.plot([toYearFraction(ydtp), toYearFraction(ydt)], \
                 [mcomps[-1], mcomps[-1]], 'g-', lw=1.5)
        
        # plt last V completeness
        #dcut = 1900
        ymax = ax.get_ylim()[-1]
        plt.plot([toYearFraction(ydt), toYearFraction(ydt)], \
                 [mcomps[-1], ymax], 'g-', lw=1.5)
        
        Npass = str(len(mvect))
        Nfail = str(len(orig_mvect) - len(mvect))
        plt.legend([h1[0], h2[0]], [Nfail+' Failed', Npass+' Passed'], loc=3, numpoints=1)
                
        tlim = ax.get_xlim()
        #dttlim = [floor(tlim[0]), ceil(tlim[1])]
        plt.xlim(tlim)
        
        # sey ylim to one
        ylims = array(ax.get_ylim())
        ylims[0] = 2.5
        plt.ylim(ylims)
        
        ###############################################################################
        # plot MFD
        ###############################################################################
        if not isnan(bval):
            # plot original data
            ax = plt.subplot(132)
            
            # plt class cum rates first        
            uidx = unique(class_cum_rates[class_idx][::-1], return_index=True, return_inverse=True)[1]
            #plt.errorbar(class_mrng[class_idx][::-1][uidx], class_cum_rates[class_idx][::-1][uidx], \
            #             yerr=[class_err_lo[class_idx][::-1][uidx], class_err_up[class_idx][::-1][uidx]], fmt='k.')
            h10 = plt.semilogy(class_mrng[class_idx][::-1][uidx], class_cum_rates[class_idx][::-1][uidx], \
                               'o', mec='b', mfc='none', ms=7)
            
            # plot best fit
            mpltmin_best = 2.0# + bin_width/2.
            plt_width = 0.1
            
            # get class betacurve
            classbetacurve, mfd_mrng = get_oq_incrementalMFD(beta, class_fn0[class_idx], \
                                                             mpltmin_best, class_mrng[class_idx][-1], \
                                                             plt_width)
            # plt class beta
            h30 = plt.semilogy(mfd_mrng, classbetacurve, '--', c='0.2')
            
            #################################################################################
            # get area normalised rates
            areanormcurve = classbetacurve * src_area[-1] / class_area[class_idx]
            h40 = plt.semilogy(mfd_mrng, areanormcurve, '--', c='r')
            
            #################################################################################
            # now plt unique values for current source
            uidx = unique(cum_rates[::-1], return_index=True, return_inverse=True)[1]
            plt.errorbar(mrng[::-1][uidx], cum_rates[::-1][uidx], \
                         yerr=[err_lo[::-1][uidx], err_up[::-1][uidx]], fmt='k.')
            h0 = plt.semilogy(mrng[::-1][uidx], cum_rates[::-1][uidx], 'ro', ms=7)
            
            # get betacurve for source
            betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mpltmin_best, mrng[-1], plt_width)
            
            h3 = plt.semilogy(mfd_mrng, betacurve, 'k-')
            
            # plot upper and lower curves
            h1 = plt.semilogy(bc_mrng_up, bc_up173, '-', c='limegreen')
            h2 = plt.semilogy(bc_mrng_up, bc_up100, 'b-')
            h4 = plt.semilogy(bc_mrng_lo, bc_lo100, 'b-')
            h5 = plt.semilogy(bc_mrng_lo, bc_lo173, '-', c='limegreen')
            
            plt.ylabel('Cumulative Rate (/yr)')
            plt.xlabel('Magnitude (MW)')
            
            # get plotting limits
            plt.xlim([2.0, bc_mrng_up[-1]+bin_width])
            yexpmin = min(floor(log10(hstack((bc_up173, bc_up100, betacurve, bc_lo100, bc_lo173)))))
            yexpmax = ceil(log10(max(class_cum_rates[class_idx])))
            #plt.ylim([10**yexpmin, 10**yexpmax])
            plt.ylim([1E-6, 100]) # for comparison across zones
            
            ###############################################################################
            # get legend text
            ###############################################################################
            
            up173_txt = '\t'.join(('Upper 1.73x', str('%0.3e' % N0_up173), str('%0.2f' % (beta-sigbeta173)), \
                                   str('%0.3f' % beta2bval(beta-sigbeta173)), str('%0.1f' % src_mmax_u[i]))).expandtabs()
                                   
            up100_txt = '\t'.join(('Upper 1.00x', str('%0.3e' % N0_up100), str('%0.2f' % (beta-sigbeta)), \
                                   str('%0.3f' % beta2bval(beta-sigbeta)), str('%0.1f' % src_mmax_u[i]))).expandtabs()
                                   
            best_txt  = '\t'.join(('Best Estimate', str('%0.3e' % fn0), str('%0.2f' % beta), \
                                   str('%0.3f' % bval), str('%0.1f' % src_mmax[i]))).expandtabs()
            
            lo100_txt = '\t'.join(('Lower 1.00x', str('%0.3e' % N0_lo100), str('%0.2f' % (beta+sigbeta)), \
                                   str('%0.3f' % beta2bval(beta+sigbeta)), str('%0.1f' % src_mmax_l[i]))).expandtabs()
                                   
            lo173_txt = '\t'.join(('Lower 1.73x', str('%0.3e' % N0_lo173), str('%0.2f' % (beta+sigbeta173)), \
                                   str('%0.3f' % beta2bval(beta+sigbeta173)), str('%0.1f' % src_mmax_l[i]))).expandtabs()
            
            # set legend title
            title = '\t'.join(('','','N Earthquakes: '+str(sum(n_obs)))).expandtabs() + '\n' \
                      + '\t'.join(('','','Regression Mmin: '+str(str(src_mmin_reg[i])))).expandtabs() + '\n\n' \
                      + '\t'.join(('','','','N0','Beta','bval','Mx')).expandtabs()
            '''
            if L08_b == True:
                title = '\t'.join(('','','N Earthquakes: '+str(sum(n_obs)))).expandtabs() + '\n' \
                          + '\t'.join(('','','Regression Mmin: '+str(str(src_mmin_reg[i])))).expandtabs() + '\n\n' \
                          + '\t'.join(('','','','N0','Beta','bval (L08)','Mx')).expandtabs()
            elif Aki_ML == True:
                title = '\t'.join(('','','N Earthquakes: '+str(sum(n_obs)))).expandtabs() + '\n' \
                          + '\t'.join(('','','Regression Mmin: '+str(str(src_mmin_reg[i])))).expandtabs() + '\n\n' \
                          + '\t'.join(('','','','N0','Beta','bval (Aki)','Mx')).expandtabs()
            else:
                title = '\t'.join(('','','N Earthquakes: '+str(sum(n_obs)))).expandtabs() + '\n' \
                          + '\t'.join(('','','Regression Mmin: '+str(str(src_mmin_reg[i])))).expandtabs() + '\n\n' \
                          + '\t'.join(('','','','N0','Beta','bval','Mx')).expandtabs()
            '''      
            leg = plt.legend([h1[0], h2[0], h3[0], h4[0], h5[0]], [up173_txt, up100_txt, best_txt, lo100_txt, lo173_txt], \
                       fontsize=9, loc=3, title=title)
            plt.setp(leg.get_title(),fontsize=9)
            
            # plt second legend
            plt.legend([h0[0], h10[0], h30[0], h40[0]], ['Zone Rates', 'Class Rates', 'Class Fit', 'Area Norm'], \
                       fontsize=10, loc=1, numpoints=1)
                       
            # replot first legend
            plt.gca().add_artist(leg)
            
            plt.grid(which='both', axis='both')
        
        ###############################################################################
        # make map
        ###############################################################################
        print 'Making map for:', src_code[i]
        ax = plt.subplot(236)
        res = 'l'
        
        '''    
        # get map bounds- for single zone
        bnds = poly.bounds
        '''
        
        # set national-scale basemap
        llcrnrlat = -44
        urcrnrlat = -6
        llcrnrlon = 107
        urcrnrlon = 152
    
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
        lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
        
        if urcrnrlat > 90.0:
            urcrnrlat = 90.0
        '''
        # get parallel/meridian spacing
        if bnds[2] - bnds[0] < 0.5:
            ll_space = 0.25
        elif bnds[2] - bnds[0] < 1.0:
            ll_space = 0.5
        elif bnds[2] - bnds[0] < 4.0:
            ll_space = 1.0
        elif bnds[2] - bnds[0] < 10.0:
            ll_space = 4.0
        else:
            ll_space = 6.0
        '''
        # draw parallels and meridians.
        ll_space = 10
        
        # set map projection
        m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                    resolution=res,area_thresh=10000.)
        
        # annotate
        m.drawcoastlines(linewidth=0.5,color='k')
        m.drawcountries()
        m.drawstates()
        m.fillcontinents(color='0.8', lake_color='1.0')
            
        # draw parallels and meridians.
        m.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
        m.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
        
        # plt earthquakes boundaries for all zones in class
        for code in class_codes[class_idx]:
            drawoneshapepoly(m, plt, sf, 'CODE', code, lw=1.5, col='b')
        
        # plt current source zone boundary
        drawoneshapepoly(m, plt, sf, 'CODE', src_code[i], lw=1.5, col='r')
        
        # map all earthquakes
        #x, y = m(eqlo, eqla)
        #size = 2.5 * eqmw
        #m.scatter(x, y, marker='o', s=size, edgecolors='0.3', mfc='none', lw=0.5)
        
        '''
        # map earthquakes that pass completeness
        for e in ev_dict:
            ms = 2.0 * e['prefmag']
            x, y = m(e['lon'], e['lat'])
            m.plot(x, y, 'o', mec='k', mfc='none', mew=1.0, ms=ms)
        '''
        # map all earthquakes in class    
        for la, lo, mag in zip(class_lavect[class_idx], class_lovect[class_idx], class_mvect[class_idx]):
            if mag >= 3.0:
                ms = 1.75 * mag
                x, y = m(lo, la)
                m.plot(x, y, 'o', mec='k', mfc='none', mew=0.5, ms=ms)
            
        # make legend - set dummy x, y params
        x, y= -100000, -100000
        
        # get marker sizes
        ms = 1.75 * 3.0
        l1 = m.plot(x, y, 'ko', ms=ms)
        ms = 1.75 * 5.0
        l2 = m.plot(x, y, 'ko', ms=ms)
        ms = 1.75 * 7.0
        l3 = m.plot(x, y, 'ko', ms=ms)
        plt.legend([l1[0], l2[0], l3[0]],['3.0', '5.0', '7.0'], fontsize=10, numpoints=1)
        
        ###############################################################################
        # make depth histogram
        ###############################################################################
        
        ax = plt.subplot(233)
        
        deprng = arange(0, 73, 4)
        
        # get depth data
        all_dep = []
        free_dep = []
        for ev in ev_dict:
            if isnan(ev['dep']) == False:
                all_dep.append(ev['dep'])
                
                if ev['fixdep'] != 1:
                    free_dep.append(ev['dep'])
        
        # first plt all data
        plt.hist(array(all_dep), bins=deprng, facecolor='w', label='Fixed Depths')
        
        # plt data with free depths
        plt.hist(array(free_dep), bins=deprng, facecolor='seagreen', label='Free Depths')
        
        # make pretty
        plt.xlabel('Hypocentral Depth (km)')
        plt.ylabel('Count')
        plt.legend()        
        
        ###############################################################################
        # make cummulative M >= 3 plot of non filtered events
        ###############################################################################
        
        ax = plt.subplot(234)
        
        # get ndays
        #td = ev_dict[-1]['datetime'] - ev_dict[0]['datetime']
        td = orig_tvect[-1] - orig_tvect[0]
        
        ndays = timedelta2days_hours_minutes(td)[0]
        
        # get events M >= 3
        dates_ge_3 = []
        dcut = 1900
        for omag, otime in zip(orig_mvect, orig_tvect):
            #if ev['prefmag'] >= 3.5 and ev['datetime'].year >= dcut:
            if omag >= 3.5 and otime.year >= dcut:
                # get decimal years
                #dates_ge_3.append(ev['datetime'].year + float(ev['datetime'].strftime('%j'))/365.) # ignore leap years for now
                dates_ge_3.append(otime.year + float(otime.strftime('%j'))/365.) # ignore leap years for now
        
        dates_ge_3 = array(dates_ge_3)        
        
        # make cummulative plot
        didx = where(dates_ge_3 > dcut)[0]
        if ndays > 0 and len(didx) > 0:
            plt.hist(dates_ge_3[didx], ndays, histtype='step', cumulative=True, color='k', lw=1.5)
            plt.xlabel('Event Year')
            plt.ylabel('Count | MW >= 3.5')
        
            # set xlims
            tlim = [int(round(x)) for x in tlim] # converting to ints
            
            plt.xlim(tlim)
            
            # sey ylim to zero
            ylims = array(ax.get_ylim())
            ylims[0] = 0
            plt.ylim(ylims)
        
        ###############################################################################
        # make src folder
        ###############################################################################
        
        srcfolder = path.join(outfolder, src_code[i])
        
        # check to see if exists
        if path.isdir(srcfolder) == False:
            mkdir(srcfolder)
        
        ###############################################################################
        # export rates file
        ###############################################################################
        # get beta curve again at consistent mags
        mpltmin_best = 2.0 + bin_width/2.
        plt_width = 0.1
        betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mpltmin_best, mrng[-1], plt_width)
        
        header = 'MAG,N_OBS,N_CUM,BIN_RTE,CUM_RTE,MFD_FIT'
        
        rate_txt = header + '\n'
        for mr in range(0,len(mrng)):
            for bm in range(0, len(mfd_mrng)):
                if around(mfd_mrng[bm], decimals=2) == around(mrng[mr], decimals=2):
                    beta_curve_val = betacurve[bm]
                    
            line = ','.join((str(mrng[mr]), str(n_obs[mr]), str(cum_num[mr]), \
                             str('%0.4e' % bin_rates[mr]), str('%0.4e' % cum_rates[mr]), \
                             str('%0.4e' % beta_curve_val))) + '\n'
            rate_txt += line
                
        # export to file
        ratefile = path.join(srcfolder, '_'.join((src_code[i], 'rates.csv')))
        f = open(ratefile, 'wb')
        f.write(rate_txt)
        f.close()
                                 
        ###############################################################################
        # export ggcat for source zone
        ###############################################################################
        catfile = path.join(srcfolder, '_'.join((src_code[i], 'passed.dat')))
        ggcat2ascii(ev_dict, catfile)
        
        # reorder out dict 
        ordidx = argsort(argsort(out_idx))
        new_dict = []
        for o in range(0, len(ordidx)):
            idx = where(ordidx == o)[0][0]
            new_dict.append(ev_out[idx])
        
        catfile = path.join(srcfolder, '_'.join((src_code[i], 'failed.dat')))
        ggcat2ascii(new_dict, catfile)
        
        ###############################################################################
        # export ggcat files to shp
        ###############################################################################
        
        # write "in" shapefile
        incat = '_'.join((src_code[i], 'passed'))+'.dat'
        inshp = '_'.join((src_code[i], 'passed'))+'.shp'
        catfile = path.join(srcfolder, incat)
        shppath = path.join(srcfolder, inshp)
        #cat2shp(catfile, shppath)
        
        # write "out" shapefile
        outcat = '_'.join((src_code[i], 'failed'))+'.dat'
        outshp = '_'.join((src_code[i], 'failed'))+'.shp'
        catfile = path.join(srcfolder, outcat)
        shppath = path.join(srcfolder, outshp)
        #cat2shp(catfile, shppath)
        
        ###############################################################################
        # save figures
        ###############################################################################
        
        # add plot sup title
        suptitle = src_name[i] + ' (' + src_code[i] + ')'
        if L08_b == True:
            suptitle += ' - L08 b-values'
        elif Aki_ML == True:
            suptitle += ' - Aki ML'
        elif Weichert == True:
            suptitle += ' - Weichert'
        else:
            suptitle += ' - Class: '+source_class
                
        plt.suptitle(suptitle, fontsize=18)
        
        pngfile = '.'.join((src_code[i], 'mfd', 'png'))
        pngpath = path.join(srcfolder, pngfile)
        plt.savefig(pngpath, format='png', bbox_inches='tight')
        
        pdffile = '.'.join((src_code[i], 'mfd', 'pdf'))
        pdfpath = path.join(srcfolder, pdffile)
        print 'Saving file:', pdfpath
        plt.savefig(pdfpath, format='pdf', bbox_inches='tight')  # causing program to crash (sometimes) on rhe-compute for unknown reason
        
        if single_src == True:
            plt.show()
        else:
            #plt.gcf().clear()
            plt.clf()
            plt.close()
    
###############################################################################
# write shapes to new shapefile
###############################################################################

# Re-read original N0 values - not sure why i need to do this,
# but seems to be getting overwritten somewhere!!!
src_n0 = get_field_data(sf, 'N0_BEST', 'float')

# get original shapefile records, and rewrite
records = sf.records()

# set shapefile to write to
w = shapefile.Writer(shapefile.POLYGON)
w.field('SRC_NAME','C','100')
w.field('CODE','C','10')
w.field('SRC_TYPE','C','10')
#w.field('CLASS','C','10')
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
w.field('MCOMP','C','30')
w.field('YMAX','F', 8, 0)
w.field('TRT','C','100')
w.field('DOMAIN','F', 2, 0)
w.field('CAT_FILE','C','50')

# make array of output field 
fields = [x[0] for x in w.fields]

# loop through original records
i = 0
for record, shape in zip(records, shapes):

    # set shape polygon
    w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
    
    # loop thru fields and match with original shapefile
    for j, field in enumerate(fields):
       
        # get field index from old shpfile
        idx = get_field_index(sf, field)
    
        # make record tuple from input shapefile
        if j == 0:
            newrec = [record[idx]]
        else:
            newrec.append(record[idx])
    
    # write new records
    # update values   
    if src_n0[i] != new_n0_b[i]:
        w.record(newrec[0], newrec[1], newrec[2], newrec[3], newrec[4], newrec[5], newrec[6], \
                 newrec[7], newrec[8], newrec[9], newrec[10], newrec[11], \
                 new_n0_b[i], new_n0_l[i], new_n0_u[i], new_bval_b[i], new_bval_l[i], new_bval_u[i], \
                 newrec[18], newrec[19], newrec[20], newrec[21], newrec[22], newrec[23], newrec[24], newrec[25])
    
    # don't edit values
    else:
        w.record(newrec[0], newrec[1], newrec[2], newrec[3], newrec[4], newrec[5], newrec[6], \
                 newrec[7], newrec[8], newrec[9], newrec[10], newrec[11], \
                 newrec[12], newrec[13], newrec[14], newrec[15], newrec[16], newrec[17], \
                 newrec[18], newrec[19], newrec[20], newrec[21], newrec[22], newrec[23], newrec[24], newrec[25])

    i += 1  
    
# now save area shapefile
newshp = path.join(rootfolder,'shapefiles',outsrcshp)
w.save(newshp)

# write projection file in WGS84
prjfile = path.join(rootfolder,'shapefiles',outsrcshp.strip().split('.shp')[0]+'.prj')
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()

###############################################################################
# map b-value
###############################################################################

# set figure
plt.clf()
fig = plt.figure(221, figsize=(13, 9))

# set national-scale basemap
llcrnrlat = -45
urcrnrlat = -5
llcrnrlon = 105
urcrnrlon = 155
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

m2 = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
            resolution='l',area_thresh=1000.)

# annotate
m2.drawcoastlines(linewidth=0.5,color='0.25')
m2.drawcountries()
m2.drawstates()
#m2.fillcontinents(color='0.8', lake_color='1.0')
    
# draw parallels and meridians.
ll_space = 6
m2.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
m2.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)

# get colour index
ncolours=12
b_min = 0.7
b_max = 1.3

cindex = []

# loop thru b values
for b in new_bval_b:
    idx = interp(b, [b_min, b_max], [0, ncolours-1])
    cindex.append(int(round(idx)))
    
# get cmap
cmap = plt.get_cmap('YlOrRd', ncolours)

# plt source zone boundary
drawshapepoly(m2, plt, sf, cindex=cindex, cmap=cmap, ncolours=ncolours, fillshape=True)

# label polygons
labelpolygon(m2, plt, sf, 'CODE')

# set colorbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.33,0.05,0.34,0.02]) # setup colorbar axes.
norm = colors.Normalize(vmin=b_min, vmax=b_max)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
#linticks = array([0.01, 0.03, 0.1, 0.3 ])
ticks = arange(b_min, b_max+0.1, 0.1)
cb.set_ticks(ticks)
labels = [str('%0.1f' % x) for x in ticks]

#cb.set_ticklabels(labels, fontsize=10)
cb.ax.set_xticklabels(labels, fontsize=10)
cb.set_label('b-value', fontsize=12)

# set filename
modelsplit = path.split(rootfolder)[-1] 
bmap = path.join(rootfolder,modelsplit+'_b_val_map.pdf')
plt.savefig(bmap, format='pdf', bbox_inches='tight')
plt.gcf().clear()
plt.clf()
plt.close()

###############################################################################
# map activity rate of M5
###############################################################################

# set figure
plt.clf()
fig = plt.figure(222, figsize=(13, 9))

# set national-scale basemap
m2 = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
            resolution='l',area_thresh=1000.)

# annotate
m2.drawcoastlines(linewidth=0.5,color='0.25')
m2.drawcountries()
m2.drawstates()
    
# draw parallels and meridians.
m2.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
m2.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)

# get M5 rates
new_beta = bval2beta(array(new_bval_b))
src_mmax = array(src_mmax)

m5_rates = array(new_n0_b) * exp(-new_beta  * 5.0) * (1 - exp(-new_beta * (src_mmax - 5.0))) \
           / (1 - exp(-new_beta * src_mmax))

# get area (in km**2) of sources for normalisation
src_area= array(src_area)
    
# normalise M5 rates by area
lognorm_m5_rates = log10(100**2 * m5_rates / src_area)
#norm_m5_rates = m5_rates
    
# get colour index
ncolours=20
r_min = -4.0
r_max = -1.5
r_rng = arange(r_min, r_max+0.1, 0.5)

cindex = []
# loop thru rates and get c-index
for r in lognorm_m5_rates:
    idx = interp(r, [r_min, r_max], [0, ncolours-1])
    cindex.append(int(round(idx)))
    
# get cmap
cmap = plt.get_cmap('rainbow', ncolours)

# plt source zone boundary
drawshapepoly(m2, plt, sf, cindex=cindex, cmap=cmap, ncolours=ncolours, fillshape=True)

# label polygons
labelpolygon(m2, plt, sf, 'CODE')

# set colorbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.33,0.05,0.34,0.02]) # setup colorbar axes.
norm = colors.Normalize(vmin=r_min, vmax=r_max)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = arange(r_min, r_max+0.5, 0.5)
cb.set_ticks(ticks)
labels = [str('%0.1e' % 10**x) for x in ticks]
cb.ax.set_xticklabels(labels, fontsize=10)
cb.set_label('M 5.0 / yr / 10,000 $\mathregular{km^{2}}$', fontsize=12)

# set filename
rmap = path.join(rootfolder,modelsplit+'_m5_rate_map.pdf')
plt.savefig(rmap, format='pdf', bbox_inches='tight')
plt.gcf().clear()
plt.clf()
plt.close()

###############################################################################
# map activity rate of M6
###############################################################################

# set figure
plt.clf()
fig = plt.figure(333, figsize=(13, 9))

# set national-scale basemap
m2 = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
            resolution='l',area_thresh=1000.)

# annotate
m2.drawcoastlines(linewidth=0.5,color='0.25')
m2.drawcountries()
m2.drawstates()
    
# draw parallels and meridians.
m2.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
m2.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)

# get M6 rates
new_beta = bval2beta(array(new_bval_b))
src_mmax = array(src_mmax)

m6_rates = array(new_n0_b) * exp(-new_beta  * 6.0) * (1 - exp(-new_beta * (src_mmax - 6.0))) \
           / (1 - exp(-new_beta * src_mmax))

# get area (in km**2) of sources for normalisation
src_area= array(src_area)
    
# normalise M6 rates by area
lognorm_m6_rates = log10(100**2 * m6_rates / src_area)
    
# get colour index
ncolours=20
r_min = -5.0
r_max = -2.5
r_rng = arange(r_min, r_max+0.1, 0.5)

cindex = []
# loop thru rates and get c-index
for r in lognorm_m6_rates:
    idx = interp(r, [r_min, r_max], [0, ncolours-1])
    cindex.append(int(round(idx)))
    
# get cmap
cmap = plt.get_cmap('rainbow', ncolours)

# plt source zone boundary
drawshapepoly(m2, plt, sf, cindex=cindex, cmap=cmap, ncolours=ncolours, fillshape=True)

# label polygons
labelpolygon(m2, plt, sf, 'CODE')

# set colorbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.33,0.05,0.34,0.02]) # setup colorbar axes.
norm = colors.Normalize(vmin=r_min, vmax=r_max)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = arange(r_min, r_max+0.5, 0.5)
cb.set_ticks(ticks)
labels = [str('%0.1e' % 10**x) for x in ticks]
cb.ax.set_xticklabels(labels, fontsize=10)
cb.set_label('M 6.0 / yr / 10,000 $\mathregular{km^{2}}$', fontsize=12)

# set filename
rmap = path.join(rootfolder,modelsplit+'_m6_rate_map.pdf')
plt.savefig(rmap, format='pdf', bbox_inches='tight')
plt.gcf().clear()
plt.clf()
plt.close()

###############################################################################
# merge all pdfs to single file
###############################################################################

from PyPDF2 import PdfFileMerger, PdfFileReader

# get input files
pdffiles = []

# make out file name
pdfbase = path.split(newshp)[-1].strip('shp')+'MERGE.pdf'
combined_pdf = path.join(rootfolder, pdfbase)

for root, dirnames, filenames in walk(rootfolder):
    #for filename in filter(filenames, '.pdf'):
    for filename in filenames:
        if filename.endswith('.pdf'):
            if filename.startswith(outsrcshp.split('.shp')[0]) == False:
                # ignore results from single src file                
                if not filename.endswith('_NSHA18_MFD.pdf'):
                    if not filename.endswith('_NSHA18_MFD.MERGE.pdf'):
                        if not filename.endswith('_NSHA18_MFD_M4.MERGE.pdf'):
                            if not filename.endswith('_NSHA18_MFD_M4_GK.MERGE.pdf'):
                                print 'Adding', filename
                                pdffiles.append(path.join(root, filename))

# now merge files
merger = PdfFileMerger()                              
for pdffile in pdffiles:                            
    merger.append(PdfFileReader(file(pdffile, 'rb')))

merger.write(combined_pdf)

###############################################################################
# write summary csv
###############################################################################

# read new shapefile
sf = shapefile.Reader(newshp)
records = sf.records()
shapes  = sf.shapes()

# make header
fields = sf.fields[1:]
simpleFields = [x[0] for x in fields]
header = ','.join(simpleFields)

csvtxt = header + '\n'
for rec in records:
    newline = ','.join(rec) + '\n'
    csvtxt += newline

csvbase = path.split(newshp)[-1].strip('shp')+'csv'
combined_csv = path.join(rootfolder, csvbase)
   
f = open(combined_csv, 'wb')
f.write(csvtxt)
f.close()
"""
###############################################################################
# make source dict for OQ input writer
###############################################################################

# set model list
model = []

# loop thru recs and make dict
for rec, shape in zip(records, shapes):
    if not float(rec[15]) == -99:
        min_mag = src_mmin[0]
        m = {'src_name':rec[0], 'src_code':rec[1], 'src_type':rec[2], 'trt':rec[23], 'src_shape':array(shape.points), 'src_dep':[float(rec[4]), float(rec[5]), float(rec[6])],
             'src_N0':[float(rec[12]), float(rec[13]), float(rec[14])], 'src_beta':[bval2beta(float(rec[15])), bval2beta(float(rec[16])), bval2beta(float(rec[17]))], 
             'max_mag':[float(rec[9]), float(rec[10]), float(rec[11])], 'min_mag':min_mag, 'src_weight':float(rec[3]), 'src_reg_wt':1}
        model.append(m)

# assume single source model for now
multimods = 'False'

# now write OQ file
oqpath = path.join(rootfolder)
write_oq_sourcefile(model, oqpath, oqpath, multimods, bestcurve)
"""