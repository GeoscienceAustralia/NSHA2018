# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 12:10:31 2017

@author: u56903
"""
from numpy import array, sqrt, where, nan, isnan, delete, hstack, diff, log10, isfinite
from tools.nsha_tools import toYearFraction, get_shapely_centroid
from shapely.geometry import Point, Polygon
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
from datetime import datetime
from catalogue_tools import weichert_algorithm, aki_maximum_likelihood, bval2beta
from oq_tools import get_oq_incrementalMFD, beta2bval
from mapping_tools import get_field_data
import shapefile
from os import path

###############################################################################
# function to parse modified HMTK catalogue with added fields
###############################################################################

def parse_hmtk_cat(hmtk_csv):
    
    print 'parsing HMTK catalogue...'

    # parse HMTK csv using modified version of HMTK parser
    parser = CsvCatalogueParser(hmtk_csv)
    hmtkcat = parser.read_file()
    
    # get number of earthquakes
    neq = len(hmtkcat.data['magnitude'])
    
    # reformat HMTK dict to one expected for code below
    cat = []
    for i in range(0, neq):
        # first make datestr
        try:
            if not isnan(hmtkcat.data['second'][i]):
                datestr = str(hmtkcat.data['eventID'][i]) \
                          + str('%2.2f' % hmtkcat.data['second'][i])
            else:
                datestr = str(hmtkcat.data['eventID'][i]) + '00.00'
                
            evdt = datetime.strptime(datestr, '%Y%m%d%H%M%S.%f')
        
        # if ID not date form, do it the hard way!
        except:
            datestr = ''.join((str(hmtkcat.data['year'][i]), str('%02d' % hmtkcat.data['month'][i]), 
                               str('%02d' % hmtkcat.data['day'][i]), str('%02d' % hmtkcat.data['hour'][i]),
                               str('%02d' % hmtkcat.data['minute'][i]), str('%0.2f' % hmtkcat.data['second'][i])))
        
            evdt = datetime.strptime(datestr, '%Y%m%d%H%M%S.%f')
            
        tdict = {'datetime':evdt, 'prefmag':hmtkcat.data['magnitude'][i], \
                 'lon':hmtkcat.data['longitude'][i], 'lat':hmtkcat.data['latitude'][i], \
                 'dep':hmtkcat.data['depth'][i], 'year':hmtkcat.data['year'][i], \
                 'month':hmtkcat.data['month'][i], 'fixdep':0, 'prefmagtype':'MW', \
                 'auth':hmtkcat.data['Agency'][i], 'mx_revML':hmtkcat.data['mx_revML'][i], \
                 'mx_origML':hmtkcat.data['mx_origML'][i], 'mx_origType':hmtkcat.data['mx_origType'][i],}
                 	
        cat.append(tdict)
    
    return cat, neq
  
###############################################################################
# function to get events within polygon
###############################################################################

def get_events_in_poly(cat, poly, depmin, depmax):
    '''
    mvect  = preferred MW
    mxvect = preferred original magnitudes
    tvect  = datetime array
    dec_tvect = decimal datetime
    ev_dict = event dictionary
    '''
    
    # set arrays
    mvect = []
    mxvect = []
    tvect = []
    dec_tvect = []
    ev_dict = []
    
    # now loop through earthquakes in cat
    for ev in cat:
        
        # check if pt in poly and compile mag and years
        pt = Point(ev['lon'], ev['lat'])
        if pt.within(poly) and ev['dep'] >= depmin and ev['dep'] <= depmax:
            mvect.append(ev['prefmag'])
            mxvect.append(ev['mx_origML']) # original catalogue mag for Mc model
            tvect.append(ev['datetime'])
            dec_tvect.append(toYearFraction(ev['datetime']))
            ev_dict.append(ev)
            
    return array(mvect), array(mxvect), array(tvect), array(dec_tvect), ev_dict

###############################################################################
# set confidence intervals from Table 1 in Weichert (1980)
###############################################################################

def get_confidence_intervals(n_obs, cum_rates):
    fup = array([1.84, 3.30, 4.64, 5.92, 7.16, 8.38, 9.58, 10.8, 12.0, 13.1, 14.3])
    dun = array([0.0, 0.173, 0.708, 1.37, 2.09, 2.84, 3.62, 4.42, 5.23, 6.06, 6.89])
    
    err_up = []
    err_lo = []
    for n in range(0, len(n_obs)):
        N = sum(n_obs[n:])
            
        if N <= 0:
            err_up.append(0)
            err_lo.append(0)
        elif N > 10:
            s = sqrt(1. / N)
            err_up.append((1. + s) * cum_rates[n] - cum_rates[n])
            err_lo.append(cum_rates[n] - (1. - s) * cum_rates[n])
        elif n < len((n_obs)):
            err_up.append(fup[N]/N * cum_rates[n] - cum_rates[n])
            err_lo.append(cum_rates[n] - dun[N]/N*cum_rates[n])
        else:
            err_up.append(0)
            err_lo.append(0)
             
    return array(err_up), array(err_lo)
    


###############################################################################
# get annualised rates
###############################################################################
def get_annualised_rates(mcomps, ycomps, mvect, mrng, bin_width, ymax):
    '''
    mcomps: list of completeness magnitudes
    ycomps: list of completeness years
    mvect:  list of event magnitudes
    mrng:   mag range of interest
    bin_width: mfd nin width
    ymax:   max year of catalogue + 1
    '''
    #print 'mcomps', mcomps
    #print 'mrng', mrng
    
    # get cumulative rates for mags >= m
    cum_num = []
    n_yrs = []
    
    # add small number to deal with precission issues
    for m in mrng:
        midx = where(mvect+1E-7 >= m)[0]
        
        # set temp cum mags & times
        cum_num.append(len(midx))
        
        # for each central mag bin, get normalisation time
        idx = where(mcomps <= m+bin_width/2)[0]
        if len(idx) > 0:
            n_yrs.append(ymax - ycomps[idx[-1]])
        else:
            n_yrs.append(0)
    
    # get events per mag bin
    n_obs = []
    for r in range(0, len(cum_num)-1):
        
        n_obs.append(cum_num[r] - cum_num[r+1])
    n_obs.append(cum_num[-1])
    
    n_obs = array(n_obs)
    n_yrs = array(n_yrs)
    bin_rates = n_obs / n_yrs
    
    #print mrng, n_obs, n_yrs, bin_rates
    
    # get cumulative rates per mag bin (/yr)
    cum_rates = []
    for m in mrng:
        midx = where( array(mrng) >= m )[0]
        cum_rates.append(sum(bin_rates[midx]))
    
    #print cum_rates
    return array(cum_rates), array(cum_num), bin_rates, n_obs, n_yrs

###############################################################################
# get events that pass completeness
###############################################################################
# assume Mc based on original catalogue magnitudes
def remove_incomplete_events(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, bin_width):
    
    # first, remove NaN magnitudes - not sure why "< mcomps[0]" fails in next block
    didx = where(isnan(mvect))[0] # use mvect here, as may be some instances where MW is not calculated
    tvect = delete(tvect, didx)
    dec_tvect = delete(dec_tvect, didx)
    mvect = delete(mvect, didx)
    mxvect = delete(mxvect, didx)
    ev_dict = delete(ev_dict, didx)
    
    # second remove all events smaller than min Mc mag minus bin_width / 2
    halfbw = bin_width / 2.
    
    didx = where(mvect < mcomps[0]-halfbw)[0]
    out_idx = didx
    tvect = delete(tvect, didx)
    dec_tvect = delete(dec_tvect, didx)
    mvect = delete(mvect, didx)
    mxvect = delete(mxvect, didx)
    ev_out = array(ev_dict)[didx]
    ev_dict = delete(ev_dict, didx)
    
    # now loop thru completeness years and mags
    for yi in range(0, len(ycomps)-1):
        # convert y to datetime
        ydt = datetime(ycomps[yi], 1, 1)
        
        # find events that meet Y and M+1 thresholds
        didx = where((tvect < ydt) & (mvect < mcomps[yi+1]-halfbw))[0]
        out_idx = hstack((out_idx, didx))
        
        # now delete events
        tvect = delete(tvect, didx)
        dec_tvect = delete(dec_tvect, didx)
        mvect = delete(mvect, didx)
        mxvect = delete(mxvect, didx)
        ev_out = hstack((ev_out, array(ev_dict)[didx]))
        ev_dict = delete(ev_dict, didx)
        
    # finally remove all events at times LT min date
    ydt = datetime(ycomps[-1], 1, 1)
    didx = where(tvect < ydt)[0]
    out_idx = hstack((out_idx, didx))
    tvect = delete(tvect, didx)
    dec_tvect = delete(dec_tvect, didx)
    mvect = delete(mvect, didx)
    mxvect = delete(mxvect, didx)
    ev_out = hstack((ev_out, array(ev_dict)[didx]))
    ev_dict = delete(ev_dict, didx)
    
    return mvect, mxvect, tvect, dec_tvect, ev_dict, out_idx, ev_out

###############################################################################
# return rates
###############################################################################

def get_mfds(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, ymax, mrng, src_mmax, \
             src_mmin_reg, src_bval_fix, src_bval_fix_sd, bin_width, poly):
    
    # remove incomplete events based on original preferred magnitudes (mxvect)
    mvect, mxvect, tvect, dec_tvect, ev_dict, out_idx, ev_out = \
         remove_incomplete_events(mvect, mxvect, tvect, dec_tvect, ev_dict, mcomps, ycomps, bin_width)
        
    # get annualised rates using preferred MW (mvect)
    cum_rates, cum_num, bin_rates, n_obs, n_yrs = \
        get_annualised_rates(mcomps, ycomps, mvect, mrng, bin_width, ymax)
            
    ###############################################################################
    # calculate MFDs if at least 50 events
    ###############################################################################
    
    # get index of min reg mag and valid mag bins
    diff_cum = abs(hstack((diff(cum_rates), 0.)))
    midx = where((mrng >= src_mmin_reg-bin_width/2.) & (isfinite(diff_cum)))[0]
    
    # check if length of midx = 0 and get highest non-zero mag
    if len(midx) == 0:
        midx = [where(isfinite(diff_cum))[0][-1]]
    
    # make sure there is at least 4 observations for b-value calculations
    if len(midx) < 5:
        idxstart = midx[0] - 1
        
        while idxstart >= 0 and len(midx) < 5:
            # if num observations greater than zero, add to midx
            if n_obs[idxstart] > 0:
                midx = hstack((idxstart, midx))
                print '    get lower mag T', midx
                
            idxstart -= 1
        
    # first, check if using fixed bval and fit curve using to solve for N0
    if src_bval_fix > 0:
        print '    Using fixed b-value =', src_bval_fix, src_bval_fix_sd        
        
        # set source beta
        bval = src_bval_fix
        beta = bval2beta(bval)
        sigb = src_bval_fix_sd
        sigbeta = bval2beta(sigb)
        
        # get dummy curve
        dummyN0 = 1.
        m_min_reg = src_mmin_reg + bin_width/2.
        bc_tmp, bc_mrng = get_oq_incrementalMFD(beta, dummyN0, m_min_reg, src_mmax, bin_width)
        
        # fit to lowest mahnitude considered
        bc_lo100 = cum_rates[midx][0] * (bc_tmp / bc_tmp[0])
        
        # scale for N0
        fn0 = 10**(log10(bc_lo100[0]) + beta2bval(beta)*bc_mrng[0])

    # do Aki ML first if N events less than 50
    elif len(mvect) >= 50 and len(mvect) < 80:
            
        # do Aki max likelihood
        bval, sigb = aki_maximum_likelihood(mrng[midx]+bin_width/2, n_obs[midx], 0.) # assume completeness taken care of
        beta = bval2beta(bval)
        sigbeta = bval2beta(sigb)
        
        # now recalc N0
        dummyN0 = 1.

        bc_tmp, bc_mrng = get_oq_incrementalMFD(beta, dummyN0, mrng[0], src_mmax, bin_width)
        
        # fit to lowest magnitude considered and observed
        Nminmag = cum_rates[midx][0] * (bc_tmp / bc_tmp[0])
        
        # !!!!!! check into why this must be done - I suspect it may be that there is an Mmax eq in the zones !!!!
        fidx = midx[0]
        
        # solve for N0
        fn0 = 10**(log10(Nminmag[0]) + bval*bc_mrng[fidx])
        
        print '    Aki ML b-value =', bval, sigb
                        
    # do Weichert for zones with more events
    elif len(mvect) >= 80:
        # calculate weichert
        bval, sigb, a_m, siga_m, fn0, stdfn0 = weichert_algorithm(array(n_yrs[midx]), \
                                               mrng[midx]+bin_width/2, n_obs[midx], mrate=0.0, \
                                               bval=1.1, itstab=1E-4, maxiter=1000)
        
        beta = bval2beta(bval)
        sigbeta = bval2beta(sigb)
    
        print '    Weichert b-value = ', bval, sigb
                    
    ###############################################################################
    # calculate MFDs using NSHA13_Background if fewer than 50 events
    ###############################################################################
    
    else:
        print 'Getting b-value from NSHA Background...'            
        # set B-value to nan
        bval = nan            
        
        # load Leonard zones
        lsf = shapefile.Reader(path.join('shapefiles','NSHA13_Background','NSHA13_Background_NSHA18_MFD.shp'))
        
        # get Leonard polygons
        l08_shapes = lsf.shapes()
        
        # get Leonard b-values
        lbval  = get_field_data(lsf, 'BVAL_BEST', 'str')
        
        # get centroid of current poly
        clon, clat = get_shapely_centroid(poly)
        point = Point(clon, clat)            
        
        # loop through zones and find point in poly
        for zone_bval, l_shape in zip(lbval, l08_shapes):
            l_poly = Polygon(l_shape.points)
            
            # check if leonard centroid in domains poly
            if point.within(l_poly):
                bval = float(zone_bval)
                
        # for those odd sites outside of L08 bounds, assign b-vale
        if isnan(bval):
            bval = 0.85
        
        beta = bval2beta(bval)
        sigb = 0.1
        sigbeta = bval2beta(sigb)
        
        # solve for N0
        fn0 = fit_a_value(bval, mrng, cum_rates, src_mmax, bin_width, midx)
        
        print '    Leonard2008 b-value =', bval, sigb
        
    # get confidence intervals        
    err_up, err_lo = get_confidence_intervals(n_obs, cum_rates)
        
    return bval, beta, sigb, sigbeta, fn0, cum_rates, ev_out, err_up, err_lo

###############################################################################
# just fit a-value
###############################################################################

def fit_a_value(bval, mrng, cum_rates, src_mmax, bin_width, midx):

     beta = bval2beta(bval)
     
     # get dummy curve
     dummyN0 = 1.
     #m_min_reg = src_mmin_reg[i] + bin_width/2.
     bc_tmp, bc_mrng = get_oq_incrementalMFD(beta, dummyN0, mrng[0], src_mmax, bin_width)
     
     # fit to lowest magnitude considered and observed
     Nminmag = cum_rates[midx][0] * (bc_tmp / bc_tmp[0])
     
     # solve for N0
     fn0 = 10**(log10(Nminmag[0]) + bval*bc_mrng[midx][0])
     
     return fn0
     
###############################################################################
# AUS ML-MW conversion
###############################################################################

# ML to MW conversion from Ghasemi (2017)
def aus_ml2mw(ml):
    a1 = 0.66199378
    a2 = 1.2156352
    a3 = 1.2156352
    mx = 4.5
    my = a1 * mx + a2
    
    if ml <= mx:
        mw = a1 * ml + a2
    else:
        mw = a3 * (ml-mx) + my
    
    return mw