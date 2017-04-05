# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 10:23:53 2017

Script to evaluate and plot ML corrections applied to the 2012 NSHA catalogue

@author: u56903
"""
from catalogue.parsers import parse_NSHA2012_catalogue
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from numpy import array, where
from os import path
import shapefile

##################################################################
# parse shapefile
##################################################################

shpfile = 'aust_boundary.shp'

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)
shape = sf.shapes()
polygon = Polygon(shape[0].points)

##################################################################
# parse shapefile
##################################################################

nsha2012csv = path.join('..', 'data', 'AUSTCAT.MW.V0.11.csv')
nsha_dict = parse_NSHA2012_catalogue(nsha2012csv)

# set arrays to fill
omag = []
rmag = []
time = []
year = []

# loop through events and extract those in Australian continent
# getting points in poly
for ev in nsha_dict:    
    # check if pt in poly and compile mag and years
    pt = Point(ev['lon'], ev['lat'])
    if pt.within(polygon):        
        omag.append(ev['mx_orig'])
        rmag.append(ev['mx_rev_ml'])
        time.append(ev['datetime'])
        year.append(ev['datetime'].year)
        
omag = array(omag)
rmag = array(rmag)
time = array(time)
year = array(year)

##################################################################
# now plot bar charts
##################################################################

ymin = 1960
ymax = 2011

#ymin = 1930
#ymax = 1970

years = range(ymin, ymax)

# loop thru mags
mrng = [3, 4, 5]

for i, m in enumerate(mrng):
    n_omag = []
    n_rmag = []

    for yr in years:
        n_omag.append(len(where((year==yr) & (omag >= m))[0]))
        n_rmag.append(len(where((year==yr) & (rmag >= m))[0]))
    
    fig = plt.figure(i+1, figsize=(14, 9))
    ax = plt.subplot(111)
    width = 0.35
    hwidth = 0.35 / 2.
    
    #ax = plt.subplots()
    plt.bar(array(years)-width, n_omag, width, color='orangered', label='Original Magnitudes')
    plt.bar(years, n_rmag, width, color='seagreen', label='Modified Magnitudes')
    if i == 0:
        plt.ylim([0, 400])
    elif i == 1:
        plt.ylim([0, 100])
        
    plt.legend()
    
    # fmt axis
    pltyears = range(ymin, ymax, 2)
    ax.set_xticks(pltyears)
    labels = [str(x) for x in pltyears]
    plt.xticks(pltyears, labels, rotation=75, ha='center')
    plt.xlim([ymin-2, ymax+1])
    plt.xlabel('Year', fontsize=18)
    plt.ylabel(' '.join(('Number of Earthquakes M',r'$\geq$',str('%0.1f' % m))), fontsize=18)
    
    plt.savefig('cmp_omag_rmag_M'+str(m)+'_'+str(ymin)+'.png',fmt='png', bbox_inches='tight')


plt.show()













