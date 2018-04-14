# -*- coding: utf-8 -*-
"""
Created on Tue May 23 13:59:45 2017

@author: u56903
"""

#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
#from hmtk.seismicity.utils import haversine
from catalogue.parsers import parse_NSHA2018_catalogue
from misc_tools import dictlist2array,timedelta2days_hours_minutes, toYearFraction, ymd2doy
from numpy import array, where, hstack, delete, arange
import matplotlib.pyplot as plt
from datetime import datetime as dt 

#hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/AUSTCAT_V0.12_hmtk_deblast.csv'
nsha_csv = 'data/NSHA18CAT.MW.V0.1.csv'

# parse HMTK csv
#parser = CsvCatalogueParser(hmtk_csv)
#nshacat = parser.read_file()

nshacat = parse_NSHA2018_catalogue(nsha_csv)
mx_orig = dictlist2array(nshacat, 'mx_orig')
mx_rev_ml = dictlist2array(nshacat, 'mx_rev_ml')
mw_pref = dictlist2array(nshacat, 'prefmag')
evdt = dictlist2array(nshacat, 'datetime')
ev_type = dictlist2array(nshacat, 'ev_type')
lat = dictlist2array(nshacat, 'lat')
lon = dictlist2array(nshacat, 'lon')

# get indexes to delete
delidx = where((ev_type=='blast') | (ev_type=='coal'))[0]
delidx = hstack((delidx, where(lat > -12)[0]))
datelim = dt(1960, 1, 1)
delidx = hstack((delidx, where(evdt < datelim)[0]))

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
minmag = 5.

fig, ax = plt.subplots(1, figsize=(10,6.25))

plt_years = arange(1960, 2019)

n_origML = []
n_corrML  = []
for py in plt_years:
   n_origML.append(len(where((mx_orig >= minmag) & (decimal_yrs >= py) & (decimal_yrs < py+1))[0]))
   n_corrML.append(len(where((mx_rev_ml >= minmag) & (decimal_yrs >= py) & (decimal_yrs < py+1))[0]))

width = 0.35  
bar1 = plt.bar(plt_years - width/2, array(n_origML), width, color='orangered')
bar2 = plt.bar(plt_years + width/2, array(n_corrML), width, color='seagreen')

ax.set_xticks(plt_years[range(0, len(plt_years)+1, 2)])
plt.xticks(rotation=45) #, ha='right')
ax.set_yticks([0, 4, 8, 12, 16, 20])   
plt.ylabel('Number of Earthquakes ML '+r'$\geq$'+' '+str(minmag))   
plt.xlabel('Years')
leg1 = ax.legend((bar1[0], bar2[0]), ('Original ML', 'Revised ML'))
leg1.get_frame().set_alpha(1.)
plt.ylim([0, 20])
plt.xlim([1958, 2018]) 
ax.yaxis.grid(ls='--')
#plt.grid(ls='--', which='y')

'''
# get average from 1960 - 1988
av_n_1960_1988 = len(where((mx_orig >= minmag) & (decimal_yrs >= 1960) & (decimal_yrs < 1988))[0]) / 28. # years
av_n_1989_2017 = len(where((mx_orig >= minmag) & (decimal_yrs >= 1989) & (decimal_yrs < 2018))[0]) / 28. # years

plt.plot([1959.65, 1987.35], [av_n_1960_1988, av_n_1960_1988], '--', c='dodgerblue', lw=2.5, label='Original ML Average Annual Number')
plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], '--', c='dodgerblue', lw=2.5)

av_n_1960_1988 = len(where((mx_rev_ml >= minmag) & (decimal_yrs >= 1960) & (decimal_yrs < 1988))[0]) / 28. # years
av_n_1989_2017 = len(where((mx_rev_ml >= minmag) & (decimal_yrs >= 1989) & (decimal_yrs < 2018))[0]) / 28. # years

plt.plot([1959.65, 1987.35], [av_n_1960_1988, av_n_1960_1988], 'k--', lw=2.5, label='Revised ML Average Annual Number')
plt.plot([1988.65, 2017.35], [av_n_1989_2017, av_n_1989_2017], 'k--', lw=2.5)
leg2 = plt.legend(loc=2)
leg2.get_frame().set_alpha(1.)
'''

plt.gca().add_artist(leg1)
plt.savefig('mag_time_bar.png', fmt='png', bbox_inches='tight', dpi=300)
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