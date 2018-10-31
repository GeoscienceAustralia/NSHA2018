# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:59:45 2017

@author: u56903
"""

#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
#from hmtk.seismicity.utils import haversine
from catalogue.parsers import parse_altmag_hmtk_catalogue
from misc_tools import dictlist2array,timedelta2days_hours_minutes, toYearFraction, ymd2doy, remove_last_cmap_colour
from gmt_tools import cpt2colormap 
from numpy import array, where, hstack, delete, arange
import matplotlib.pyplot as plt
from datetime import datetime as dt 
import matplotlib as mpl
mpl.style.use('classic')

#hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/AUSTCAT_V0.12_hmtk_deblast.csv'
'''
nsha_csv = 'data/NSHA18CAT.MW.V0.1.csv'
nshacat = parse_NSHA2018_catalogue(nsha_csv)
'''

# get colours
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//Paired_10.cpt'
ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

# parse HMTK csv - use declustered catalogue
hmtk_csv = 'data//NSHA18CAT_V0.2_hmtk_declustered.csv'
nshacat = parse_altmag_hmtk_catalogue(hmtk_csv)[0]

'''
tmpdict = {'datetime':ev_date, 'lon':lon, 'lat':lat, 'dep':dep,
                   'mx_origML':mx_origML, 'mx_revML': mx_revML,
                   'mw_pref':mw_pref, 'mw_qds':mw_qds, 'mw_qde':mw_qde, 
                   'mw_ble':mw_ble, 'mw_alt_ble':mw_alt_ble, 
                   'mw_alt_qde':mw_alt_qde, 'prefmag':mw_pref}
'''

mx_orig = dictlist2array(nshacat, 'mx_origML')
mx_rev_ml = dictlist2array(nshacat, 'mx_revML')
mw_pref = dictlist2array(nshacat, 'prefmag')
evdt = dictlist2array(nshacat, 'datetime')
#ev_type = dictlist2array(nshacat, 'ev_type')
lat = dictlist2array(nshacat, 'lat')
lon = dictlist2array(nshacat, 'lon')

# get indexes to delete
#delidx = where((ev_type=='blast') | (ev_type=='coal'))[0]
#delidx = hstack((delidx, where(lat > -12)[0]))
datelim = dt(1960, 1, 1)
#delidx = hstack((delidx, where(evdt < datelim)[0]))
delidx = where(evdt < datelim)[0]

# delete events
mx_orig = delete(mx_orig, delidx)
mx_rev_ml = delete(mx_rev_ml, delidx)
mw_pref = delete(mw_pref, delidx)
evdt = delete(evdt, delidx)

# get decimal years
decimal_yrs = []
for ed in evdt:
    #decimal_yrs.append(toYearFraction(ed))
    decimal_yrs.append(ed.year + float(ymd2doy(ed.year, ed.month, ed.day)) / 365.)
decimal_yrs = array(decimal_yrs)    

# get plotting indices
minmags = [4.5, 5.0]
pltmean = [False, True]
pltmean = [True]

# loop thru years
i = 0
for minmag in minmags:
    fig, ax = plt.subplots(1, figsize=(14, 8))
    
    width = 0.35  
        
    # first plot filled box between 1986 & 1992
    if minmag == 5.0:
        fy =  7
    else:
        fy =  13
    fx1 = 1986 - width
    fx2 = 1992 + width
    plt.fill([fx1, fx2, fx2, fx1, fx1], [0, 0, fy, fy, 0], '0.8', edgecolor='0.8')
    plt.text(1989, 8.1, '  Australian $\mathregular{M_L}$\n'+'equations developed', rotation=90., ha='center', va='bottom', fontsize=13)
        
    
    for pm in pltmean:
        i += 1
        
        plt_years = arange(1960, 2019)
        
        n_origML = []
        n_corrML  = []
        for py in plt_years:
           n_origML.append(len(where((mx_orig >= minmag) & (decimal_yrs >= py) & (decimal_yrs < py+1))[0]))
           n_corrML.append(len(where((mx_rev_ml >= minmag) & (decimal_yrs >= py) & (decimal_yrs < py+1))[0]))
        
        #bar1 = plt.bar(plt_years - width/2, array(n_origML), width, color='orangered')
        #bar2 = plt.bar(plt_years + width/2, array(n_corrML), width, color='seagreen')
        #bar1 = plt.bar(plt_years - width/2, array(n_origML), width, color=cs[1])
        #bar2 = plt.bar(plt_years + width/2, array(n_corrML), width, color=cs[7])
        bar1 = plt.bar(plt_years - width, array(n_origML), width, color=cs[1])
        bar2 = plt.bar(plt_years + 0., array(n_corrML), width, color=cs[7])
        
        #ax.set_xticks(plt_years[range(0, len(plt_years)+1, 2)])
        xticks = plt_years[range(0, len(plt_years)+1, 2)]
        xtick_labels = [str(x) for x in xticks]
        plt.xticks(xticks)
        ax.set_xticklabels(xtick_labels)
        
        plt.xticks(rotation=65) #, ha='right')
        #ax.set_yticks([0, 4, 8, 12, 16, 20])   
        plt.ylabel('Number of Earthquakes $\mathregular{M_X}$ '+r'$\geq$'+' '+str(minmag), fontsize=17)   
        plt.xlabel('Year', fontsize=17)
        
        leg1 = ax.legend((bar1[0], bar2[0]), ('Original $\mathregular{M_{LH}}$', 'Revised $\mathregular{M_{LR}}$'))
        leg1.get_frame().set_alpha(1.)
        plt.xlim([1958, 2018]) 
        ax.yaxis.grid(ls=':')
        #plt.grid(ls='--', which='y')
        
        if minmag == 5.0:
            ax.set_yticks([0, 2, 4, 6])  
            plt.ylim([0, 7])
             
        
        else:
            
            ax.set_yticks([0, 4, 8, 12])   
            plt.ylim([0, 13])
        
        # get average from 1960 - 1988
        if pm == True:
            av_n_1960_1988 = len(where((mx_orig >= minmag) & (decimal_yrs >= 1960) & (decimal_yrs <= 1988))[0]) / 29. # years
            av_n_1989_2017 = len(where((mx_orig >= minmag) & (decimal_yrs >= 1989) & (decimal_yrs < 2018))[0]) / 28. # years
            
            plt.plot([1959.65, 1988.35], [av_n_1960_1988, av_n_1960_1988], '--', c='darkblue', lw=2.5, label='Original $\mathregular{M_{LH}}$ Average Annual Number')
            plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], '--', c='darkblue', lw=2.5)
            #plt.plot([1959.65, 1988.35], [av_n_1960_1988, av_n_1960_1988], '--', c=cs[-1], lw=2.5, label='Original ML Average Annual Number')
            #plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], '--', c=cs[-1], lw=2.5)
            
            av_n_1960_1988 = len(where((mx_rev_ml >= minmag) & (decimal_yrs >= 1960) & (decimal_yrs <= 1988))[0]) / 29. # years
            av_n_1989_2017 = len(where((mx_rev_ml >= minmag) & (decimal_yrs >= 1989) & (decimal_yrs < 2018))[0]) / 28. # years
            
            plt.plot([1959.65, 1988.35], [av_n_1960_1988, av_n_1960_1988], '--', c='orangered', lw=2.5, label='Revised $\mathregular{M_{LR}}$ Average Annual Number')
            plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], '--', c='orangered', lw=2.5)
            #plt.plot([1959.65, 1988.35], [av_n_1960_1988, av_n_1960_1988], '--', c=cs[5], lw=2.5, label='Revised ML Average Annual Number')
            #plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], '--', c=cs[5], lw=2.5)
            
            leg2 = plt.legend(loc=2)
            leg2.get_frame().set_alpha(1.)
            
        
        plt.gca().add_artist(leg1)
        
        if pm == True:
            figname = '_'.join(('mag_time_bar',str(minmag),'avplt.png'))
        else:
            figname = '_'.join(('mag_time_bar',str(minmag)+'.png'))
        plt.savefig(figname, fmt='png', bbox_inches='tight', dpi=300)
        
plt.show()


"""
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
            
"""            