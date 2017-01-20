from numpy import array, arange, argsort, sort, where, delete, hstack, sqrt, \
                  unique, mean, percentile, log10, zeros_like, ceil, floor, \
                  ones_like, isnan
from os import path, sep, mkdir, getcwd, system, walk, system
from shapely.geometry import Point, Polygon
from datetime import datetime
from sys import argv
import shapefile
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# import non-standard functions
try:
    from catalogue_tools import weichert_algorithm, bval2beta
    from fault_tools import get_oq_incrementalMFD, beta2bval, bval2beta
    from mapping_tools import get_field_data, get_field_index, drawoneshapepoly
    from catalogue.parsers import parse_ggcat
    from misc_tools import listdir_extension
    from make_nsha_oq_inputs import write_oq_sourcefile
except:
    cwd = getcwd().split(sep)
    pythonpath = sep.join(pt[0:-3])+sep+'tools'
    print '\nSet environmental variables, e.g.:\n\nexport PYTHONPATH='+pythonpath+':$PYTHONPATH\n'

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

###############################################################################
# get parameter values
###############################################################################

# load param file
lines = open(paramfile).readlines()
rootfolder  = lines[0].split('=')[-1].strip()
ggcatfile   = lines[1].split('=')[-1].strip()
dec_flag    = lines[2].split('=')[-1].strip() # decluster flag
shpfile     = lines[3].split('=')[-1].strip()
outfolder   = path.join(rootfolder, lines[4].split('=')[-1].strip())
outsrcshp   = lines[5].split('=')[-1].strip()
bin_width   = float(lines[6].split('=')[-1].strip())

# get export folder
now = datetime.now()
shpname = path.split(shpfile)
if dec_flag == 'True':
    #outfolder = '_'.join((shpname[-1][0:-4], 'DEC', now.strftime('%Y-%m-%d')))
    outfolder += '_'.join(('DEC', now.strftime('%Y-%m-%d')))
else:
    outfolder += '_' + now.strftime('%Y-%m-%d')

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
src_sheef = get_field_data(sf, 'CAT_FILE', 'str')
sortind = argsort(src_code)

# initiate new arrays for writing new shpfile
new_bval_b = src_bval  
new_bval_l = src_bval_l
new_bval_u = src_bval_u
new_n0_b = src_n0
new_n0_l = src_n0_l
new_n0_u = src_n0_u

print src_n0, new_n0_b
#new_polys = []

'''
for i, si in enumerate(sortind):
    print i, si, src_code[sortind[i]]
'''

# if single source remove unnecessary data
if single_src == True:
    srcidx = where(array(src_code) == single_zone)[0]
    ssi = where(sortind == srcidx)[0]

# do all sources except Beaufort-Mackenzie Convergence (BMC)
else:
    srcidx = range(len(src_code))
    
    
###############################################################################
# parse GGCat
###############################################################################

ggcat = parse_ggcat(ggcatfile)

###############################################################################
# loop thru zones 
###############################################################################

# loop thru source zones
for i in srcidx:

    ###############################################################################
    # get earthquakes that pass completeness
    ###############################################################################
    
    # get polygon of interest
    poly = polygons[i]
    
    print '\nCalculating MFD for:', src_code[i]
    
    mvect = []
    tvect = []
    ev_dict = []
    out_idx = []
    
    """
    # set SHEEF for given depth range
    if src_sheef[i] == 'FULL':
        ggcat = fullsheef
    elif src_sheef[i] == 'DEEP':
        ggcat = deepsheef
    elif src_sheef[i] == 'CRUST':
        ggcat = crustsheef
    """
    
    # get max decimal year and round up!
    lr = ggcat[-1]
    max_comp_yr = lr['year']+lr['month']/12.
    
    # now loop through earthquakes in SHEEF
    for s in ggcat:
        
        # check if pt in poly and compile mag and years
        pt = Point(s['lon'], s['lat'])
        if pt.within(poly):
            mvect.append(s['prefmag'])
            tvect.append(s['datetime'])
            ev_dict.append(s)
            
    # get annual occurrence rates for each mag bin
    ycomps = array([int(x) for x in src_ycomp[i].split(';')])
    mcomps = array([float(x) for x in src_mcomp[i].split(';')])
    
    # set mag bins
    mrng = arange(min(mcomps)-bin_width/2, src_mmax[i], bin_width)
    #mrng = arange(min(mcomps), src_mmax[i], bin_width)
    
    # convert lists to arrays
    mvect = array(mvect)
    tvect = array(tvect)
    
    # keep original vectors for plotting
    orig_mvect = mvect
    orig_tvect = tvect
    
    # first remove all events smaller than min Mc mag
    didx = where(mvect < mcomps[0])[0]
    out_idx = didx
    tvect = delete(tvect, didx)
    mvect = delete(mvect, didx)
    ev_out = array(ev_dict)[didx]
    ev_dict = delete(ev_dict, didx)
    
    # now loop thru completeness years and mags
    for yi in range(0, len(ycomps)-1):
        # convert y to datetime
        ydt = datetime(ycomps[yi], 1, 1)
        
        # find events that meet Y and M+1 thresholds
        didx = where((tvect < ydt) & (mvect < mcomps[yi+1]))[0]
        out_idx = hstack((out_idx, didx))
        
        # now delete events
        tvect = delete(tvect, didx)
        mvect = delete(mvect, didx)
        ev_out = hstack((ev_out, array(ev_dict)[didx]))
        ev_dict = delete(ev_dict, didx)
        
    # finally remove all events at times LT min date
    ydt = datetime(ycomps[-1], 1, 1)
    didx = where(tvect < ydt)[0]
    out_idx = hstack((out_idx, didx))
    tvect = delete(tvect, didx)
    mvect = delete(mvect, didx)  
    ev_out = hstack((ev_out, array(ev_dict)[didx]))
    ev_dict = delete(ev_dict, didx)  
    
    # skip zone if no events pass completeness
    if len(mvect) != 0:
    
    ###############################################################################
    # get annualised rates
    ###############################################################################
        
        # get cumulative rates for mags >= m
        cum_num = []
        n_yrs = []
        
        for m in mrng:
            midx = where(array(mvect) >= m)[0]
            
            # set temp cum mags & times
            cum_num.append(len(midx))
            
            
            # for each central mag bin, get normalisation time
            src_ymax[i] = max_comp_yr
            idx = where(mcomps <= m+bin_width/2)[0]
            if len(idx) > 0:
                n_yrs.append(src_ymax[i] - ycomps[idx[-1]])
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
        
        # get cumulative rates per mag bin (/yr)
        cum_rates = []
        for m in mrng:
            midx = where( array(mrng) >= m )[0]
            cum_rates.append(sum(bin_rates[midx]))
        
        cum_rates = array(cum_rates)
        
        ###############################################################################
        # calculate MFDs`
        ###############################################################################
        
        # get magnitude indices being considered for regression
        
        midx = where(mrng >= src_mmin_reg[i])[0]
        
        # if beta not fixed, do Weichert
        if src_bval_fix[i] == -99:
            
            # calculate weichert
            bval, sigb, a_m, siga_m, fn0, stdfn0 = weichert_algorithm(array(n_yrs[midx]), \
                                                   mrng[midx]+bin_width/2, n_obs[midx], mrate=0.0, \
                                                   bval=1.0, itstab=1E-5, maxiter=1000)
            
            beta = bval2beta(bval)
            sigbeta = bval2beta(sigb)
        
        # else, fit curve using fixed beta and solve for N0
        else:
            
            # set source beta
            beta = src_bval_fix[i]
            bval = beta2bval(beta)
            sigbeta = src_bval_fix_sd[i]
            
            # get dummy curve
            dummyN0 = 1.
            bc_tmp, bc_mrng = get_oq_incrementalMFD(beta, dummyN0, src_mmin_reg[i], src_mmax[i], bin_width)
            
            # fit to lowest mahnitude considered
            bc_lo100 = cum_rates[midx][0] * (bc_tmp / bc_tmp[0])
            
            # scale for N0
            fn0 = 10**(log10(bc_lo100[0]) + beta2bval(beta)*bc_mrng[0])
        
        print 'beta = ', bval2beta(bval)
        
        ###############################################################################
        # set confidence intervals from Table 1 in Weichert (1980)
        ###############################################################################
        
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
                 
                
        err_up = array(err_up)
        err_lo = array(err_lo)
        
        ###############################################################################
        # get upper and lower MFD bounds
        ###############################################################################
        sigbeta173 = 1.73 * sigbeta
        sigb173 = 1.73 * sigb
        
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
        
        print 'Filling new values for', src_code[i], '\n'
        new_bval_b[i] = bval
        new_bval_l[i] = bval-sigb173
        new_bval_u[i] = bval+sigb173
        new_n0_b[i]   = fn0
        new_n0_l[i]   = N0_lo173
        new_n0_u[i]   = N0_up173
        
        ###############################################################################
        # plot earthquakes that pass completeness
        ###############################################################################
    
        fig = plt.figure(i, figsize=(20, 9))
        
        # plot original data
        ax = plt.subplot(231)
        
        # cannot plot prior to 1700 - not sure why
        dcut = datetime(1700,1,1,0,0)
        didx = where(orig_tvect > dcut)[0]
        h1 = plt.plot(orig_tvect[didx], orig_mvect[didx], 'bo')
        
        # now loop thru completeness years and mags
        for yi in range(0, len(ycomps)-1):
            # convert y to datetime
            ydt = datetime(ycomps[yi], 1, 1)
            
            # plt H completeness ranges
            if yi == 0:
                plt.plot([ydt, max(tvect)], [mcomps[yi], mcomps[yi]], 'g-', lw=1.5)
            else:
                ydtp = datetime(ycomps[yi-1], 1, 1)
                plt.plot([ydtp, ydt], [mcomps[yi], mcomps[yi]], 'g-', lw=1.5)
            
            # plt V completeness ranges
            plt.plot([ydt, ydt], [mcomps[yi], mcomps[yi+1]], 'g-', lw=1.5)
        
        # plt last H completeness range
        ydt = datetime(ycomps[-1], 1, 1)
        if ydt < dcut:
            ydt = dcut
        ydtp = datetime(ycomps[-2], 1, 1)
        plt.plot([ydtp, ydt], [mcomps[-1], mcomps[-1]], 'g-', lw=1.5)
        
        # plt last V completeness
        ymax = ax.get_ylim()[-1]
        plt.plot([ydt, ydt], [mcomps[-1], ymax], 'g-', lw=1.5)
        
        didx = where(tvect > dcut)[0]        
        h2 = plt.plot(tvect[didx], mvect[didx], 'ro')
        plt.ylabel('Magnitude (MW)')
        plt.xlabel('Date')
        #plt.title(src_code[i] + ' Catalogue Completeness')
        
        Npass = str(len(mvect))
        Nfail = str(len(orig_mvect) - len(mvect))
        plt.legend([h1[0], h2[0]], [Nfail+' Failed', Npass+' Passed'], loc=3, numpoints=1)
        
        # get xlims
        tlim = ax.get_xlim()
        
        ###############################################################################
        # plot MFD
        ###############################################################################
        
        # plot original data
        ax = plt.subplot(132)
        
        # plt unique values
        uidx = unique(cum_rates[::-1], return_index=True, return_inverse=True)[1]
        plt.errorbar(mrng[::-1][uidx], cum_rates[::-1][uidx], \
                     yerr=[err_lo[::-1][uidx], err_up[::-1][uidx]], fmt='k.')
        h1 = plt.semilogy(mrng[::-1][uidx], cum_rates[::-1][uidx], 'ro', ms=8)
        
        # plot best fit
        mpltmin_best = 2.0
        plt_width = 0.1
        betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mpltmin_best, mrng[-1], plt_width)
        
        #betacurve, mfd_mrng = get_oq_incrementalMFD(beta, fn0, mmin, mrng[-1], bin_width)
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
        yexpmax = max(ceil(log10(betacurve)))
        plt.ylim([10**yexpmin, 10**yexpmax])
        
        ###############################################################################
        # get legend text
        ###############################################################################
        
        up173_txt = '\t'.join(('Upper 1.73x', str('%0.1f' % N0_up173), str('%0.2f' % (beta-sigbeta173)), \
                               str('%0.3f' % beta2bval(beta-sigbeta173)), str('%0.1f' % src_mmax_u[i]))).expandtabs()
                               
        up100_txt = '\t'.join(('Upper 1.00x', str('%0.1f' % N0_up100), str('%0.2f' % (beta-sigbeta)), \
                               str('%0.3f' % beta2bval(beta-sigbeta)), str('%0.1f' % src_mmax_u[i]))).expandtabs()
                               
        best_txt  = '\t'.join(('Best Estimate', str('%0.1f' % fn0), str('%0.2f' % beta), \
                               str('%0.3f' % bval), str('%0.1f' % src_mmax[i]))).expandtabs()
        
        lo100_txt = '\t'.join(('Lower 1.00x', str('%0.1f' % N0_lo100), str('%0.2f' % (beta+sigbeta)), \
                               str('%0.3f' % beta2bval(beta+sigbeta)), str('%0.1f' % src_mmax_l[i]))).expandtabs()
                               
        lo173_txt = '\t'.join(('Lower 1.73x', str('%0.1f' % N0_lo173), str('%0.2f' % (beta+sigbeta173)), \
                               str('%0.3f' % beta2bval(beta+sigbeta173)), str('%0.1f' % src_mmax_l[i]))).expandtabs()
        
        # set legend title
        title = '\t'.join(('','','N Earthquakes: '+str(sum(n_obs)))).expandtabs() + '\n' \
              + '\t'.join(('','','Regression Mmin: '+str(str(src_mmin_reg[i])))).expandtabs() + '\n\n' \
              + '\t'.join(('','','N0','Beta','bval','Mx')).expandtabs()
              
        leg = plt.legend([h1[0], h2[0], h3[0], h4[0], h5[0]], [up173_txt, up100_txt, best_txt, lo100_txt, lo173_txt], \
                   fontsize=9, loc=3, title=title)
        plt.setp(leg.get_title(),fontsize=9)
        plt.grid(which='both', axis='both')
        
        ###############################################################################
        # make map
        ###############################################################################
        
        ax = plt.subplot(236)
        res = 'i'
        
        # get map bounds
        bnds = poly.bounds
        mbuff = 0.2 * (bnds[2] - bnds[0])
        llcrnrlat = bnds[1] - mbuff/2.
        urcrnrlat = bnds[3] + mbuff/2.
        llcrnrlon = bnds[0] - mbuff
        urcrnrlon = bnds[2] + mbuff
        lon_0 = mean([llcrnrlon, urcrnrlon])
        lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
        lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
        
        if urcrnrlat > 90.0:
            urcrnrlat = 90.0
        
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
        
        # set map projection
        m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                    projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                    resolution=res,area_thresh=1000.)
        
        # annotate
        m.drawcoastlines(linewidth=0.5,color='k')
        m.drawcountries()
        m.drawstates()
        m.fillcontinents(color='0.8', lake_color='1.0')
            
        # draw parallels and meridians.
        m.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
        m.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
        
        # plt source zone boundary
        drawoneshapepoly(m, plt, sf, 'CODE', src_code[i], lw=2., col='r')
        
        # map all earthquakes
        #x, y = m(eqlo, eqla)
        #size = 2.5 * eqmw
        #m.scatter(x, y, marker='o', s=size, edgecolors='0.3', mfc='none', lw=0.5)
        
        # map earthquakes that pass completeness
        for e in ev_dict:
            ms = 2.5 * e['prefmag']
            x, y = m(e['lon'], e['lat'])
            m.plot(x, y, 'o', mec='k', mfc='none', mew=1.0, ms=ms)
            
        # make legend - set dummy x, y params
        x, y= -100000, -100000
        
        # get marker sizes
        ms = 2.5 * 3.0
        l1 = m.plot(x, y, 'ko', ms=ms)
        ms = 2.5 * 5.0
        l2 = m.plot(x, y, 'ko', ms=ms)
        ms = 2.5 * 7.0
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
        # make cummulative M >= 3 plot
        ###############################################################################
        
        ax = plt.subplot(234)
        
        # get ndays
        td = ev_dict[-1]['datetime'] - ev_dict[0]['datetime']
        
        def timedelta2days_hours_minutes(td):
            return td.days, td.seconds//3600, (td.seconds//60)%60
        
        ndays = timedelta2days_hours_minutes(td)[0]
        
        # get events M >= 3
        dcut = 1900. # date cutoff for plotting - datetime cannot handle dates earlier than 1900???
        dates_ge_3 = []
        for ev in ev_dict:
            if ev['prefmag'] >= 3.0 and ev['datetime'].year >= dcut:
                # get decimal years
                dates_ge_3.append(ev['datetime'].year + float(ev['datetime'].strftime('%j'))/365.) # ignore leap years for now
        
        dates_ge_3 = array(dates_ge_3)        
        
        # make cummulative plot
        didx = where(dates_ge_3 > dcut)[0]
        if ndays > 0:
            plt.hist(dates_ge_3[didx], ndays, histtype='step', cumulative=True, color='k', lw=1.5)
            plt.xlabel('Event Year')
            plt.ylabel('Count | MW >= 3.0')
        
            # set xlims
            tlim = [int(round(x)) for x in tlim] # converting to ints
            dttlim = [datetime.fromordinal(tlim[0]), datetime.fromordinal(tlim[1])] # converting from ordianl back to datetime
            # now convert to decimal year - very ugly!
            dttlim = [floor(dttlim[0].year), ceil(dttlim[1].year)]
            
            plt.xlim(dttlim)
        
        ###############################################################################
        # make src folder
        ###############################################################################
        
        srcfolder = path.join(outfolder, src_code[i])
        
        # check to see if exists
        if path.isdir(srcfolder) == False:
            mkdir(srcfolder)
        
        ###############################################################################
        # export ggcat for source zone
        ###############################################################################
        sheefpath = path.join(srcfolder, '_'.join((src_code[i], 'SHEEF_passed.dat')))
        #write_sheef(ev_dict, sheefpath)
        
        # reorder out dict 
        ordidx = argsort(argsort(out_idx))
        new_dict = []
        for o in range(0, len(ordidx)):
            idx = where(ordidx == o)[0][0]
            new_dict.append(ev_out[idx])
        
        sheefpath = path.join(srcfolder, '_'.join((src_code[i], 'SHEEF_failed.dat')))
        #write_sheef(new_dict, sheefpath)
        
        ###############################################################################
        # export ggcat files to shp
        ###############################################################################
        
        # write "in" shapefile
        insheef = '_'.join((src_code[i], 'SHEEF_passed'))+'.dat'
        inshp = '_'.join((src_code[i], 'SHEEF_passed'))+'.shp'
        sheefpath = path.join(srcfolder, insheef)
        shppath = path.join(srcfolder, inshp)
        #sheef2shp(sheefpath, shppath)
        
        # write "out" shapefile
        outsheef = '_'.join((src_code[i], 'SHEEF_failed'))+'.dat'
        outshp = '_'.join((src_code[i], 'SHEEF_failed'))+'.shp'
        sheefpath = path.join(srcfolder, outsheef)
        shppath = path.join(srcfolder, outshp)
        #sheef2shp(sheefpath, shppath)
        
        ###############################################################################
        # save figures
        ###############################################################################
        
        # add plot sup title
        plt.suptitle(src_name[i] + ' (' + src_code[i] + ')', fontsize=18)
        
        pngfile = '.'.join((src_code[i], 'mfd', 'png'))
        pngpath = path.join(srcfolder, pngfile)
        #plt.savefig(pngpath, format='png', bbox_inches='tight')
        
        pdffile = '.'.join((src_code[i], 'mfd', 'pdf'))
        pdfpath = path.join(srcfolder, pdffile)
        plt.savefig(pdfpath, format='pdf', bbox_inches='tight')
        
        if single_src == True:
            plt.show()
        else:
            plt.clf()
    
###############################################################################
# write shapes to new shapefile
###############################################################################

# get original shapefile records, and rewrite
records = sf.records()

# set shapefile to write to
w = shapefile.Writer(shapefile.POLYGON)
w.field('SRC_NAME','C','100')
w.field('CODE','C','10')
w.field('SRC_TYPE','C','10')
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
            #newrec = (record[idx],)
            newrec = [record[idx]]
        else:
            newrec.append(record[idx])
    
    # write new records
    print 'N0', src_n0[i], new_n0_b[i]
    if src_n0[i] != new_n0_b[i]:
        w.record(newrec[0], newrec[1], newrec[2], newrec[3], newrec[4], newrec[5], newrec[6], \
                 newrec[7], newrec[8], newrec[9], newrec[10], newrec[11], \
                 new_n0_b[i], new_n0_l[i], new_n0_u[i], new_bval_b[i], new_bval_l[i], new_bval_u[i], \
                 newrec[18], newrec[19], newrec[20], newrec[21], newrec[22], newrec[23], newrec[24], newrec[25])
    else:
        w.record(newrec[0], newrec[1], newrec[2], newrec[3], newrec[4], newrec[5], newrec[6], \
                 newrec[7], newrec[8], newrec[9], newrec[10], newrec[11], \
                 newrec[12], newrec[13], newrec[14], newrec[15], newrec[16], newrec[17], \
                 newrec[18], newrec[19], newrec[20], newrec[21], newrec[22], newrec[23], newrec[24], newrec[25])

    i += 1  
    
# now save area shapefile
newshp = path.join(outfolder,outsrcshp)
w.save(newshp)

# write projection file
prjfile = path.join(outfolder, outsrcshp.strip().split('.shp')[0]+'.prj')
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()

# not sure why, but delete "failed files from working folder
system('rm *failed*')
    
###############################################################################
# merge all pdfs to single file
###############################################################################

from PyPDF2 import PdfFileMerger, PdfFileReader 

# get input files
pdffiles = []
for root, dirnames, filenames in walk(outfolder):
    #for filename in filter(filenames, '.pdf'):
    for filename in filenames:
        if filename.endswith('.pdf'):
            print filename
            pdffiles.append(path.join(root, filename))

# now merge files
merger = PdfFileMerger()                              
for pdffile in pdffiles:                            
    merger.append(PdfFileReader(file(pdffile, 'rb')))

merger.write(newshp.strip().split('.shp')[0]+'.pdf')

###############################################################################
# make source dict for OQ input writer
###############################################################################

# read new shapefile
sf = shapefile.Reader(newshp)
records = sf.records()
shapes  = sf.shapes()

model = []

# loop thru recs and make dict
for rec, shape in zip(records, shapes):
    m = {'src_name':rec[0], 'src_code':rec[1], 'src_type':rec[2], 'trt':rec[23], 'src_shape':array(shape.points), 'src_dep':[float(rec[4]), float(rec[5]), float(rec[6])],
         'src_N0':[float(rec[12]), float(rec[13]), float(rec[14])], 'src_beta':[bval2beta(float(rec[15])), bval2beta(float(rec[16])), bval2beta(float(rec[17]))], 
         'max_mag':[float(rec[9]), float(rec[10]), float(rec[11])], 'min_mag':float(rec[7]), 'src_weight':float(rec[3]), 'src_reg_wt':1}
    model.append(m)

# assume single source model for now
multimods = 'False'

# now write OQ file
oqpath = path.join(outfolder)
write_oq_sourcefile(model, oqpath, oqpath, multimods)