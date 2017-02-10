# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 13:30:27 2017

@author: u56903
"""

from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser#, CsvCatalogueWriter
from hmtk.seismicity.utils import decimal_year, haversine
from parsers import parse_ggcat
from writers import ggcat2hmtk_csv
import numpy as np
import datetime as dt

from os import path

ggcatcsv = path.join('data', 'GGcat-161025.csv')
ggdict = parse_ggcat(ggcatcsv)

# set HMTK file name
hmtk_csv = ggcatcsv.split('.')[0] + '_hmtk.csv'

# write HMTK csv
ggcat2hmtk_csv(ggdict, hmtk_csv)

# parse HMTK csv
parser = CsvCatalogueParser(hmtk_csv)
ggcat = parser.read_file()


# try following HMTK format
print  'Using Leonard 2008 declustering...'
catalogue = ggcat

# Get relevant parameters
neq = len(catalogue.data['magnitude'])  # Number of earthquakes
# Get decimal year (needed for time windows)
year_dec = decimal_year(
    catalogue.data['year'], catalogue.data['month'],
    catalogue.data['day'])
# Get space and time windows corresponding to each event
'''
sw_space, sw_time = (
    config['time_distance_window'].calc(catalogue.data['magnitude']))

# Initial Position Identifier
eqid = np.arange(0, neq, 1)
# Pre-allocate cluster index vectors
vcl = np.zeros(neq, dtype=int)
# Sort magnitudes into descending order
id0 = np.flipud(np.argsort(catalogue.data['magnitude'],
                           kind='heapsort'))
longitude = catalogue.data['longitude'][id0]
latitude = catalogue.data['latitude'][id0]
#sw_space = sw_space[id0]
#sw_time = sw_time[id0]
year_dec = year_dec[id0]
eqid = eqid[id0]
'''

method = 'Leonard08'
# set periods of confidence
test_day_1 = dt.datetime(1960,1,1) # Earthquakes older than this are assumed to be very poorly located.
test_day_2 = dt.datetime(1970,1,1) # Earthquakes older than this are assumed to be poorly located.
delta_mag = 0.989 # Aftershocks must be less than 0.989 m.u. of the main shock (95% of the moment).

# set time window
if method == 'Leonard08':
    max_time = 10**((catalogue.data['magnitude']-1.85)*0.69)
elif method == 'Stien08':
    max_time = 10**((catalogue.data['magnitude']-2.7)*1.1) + 4.0
    
# get event time datevector
evdate = []
for i in range(0, neq):

    # get event dattime
    evdate.append(dt.datetime(catalogue.data['year'][i], \
                              catalogue.data['month'][i], \
                              catalogue.data['day'][i]))

evdate = np.array(evdate)

# set flag for events to delete
flagvector = np.zeros(neq, dtype=int)

# loop through earthquakes
print 'Looping through events...'
for i in range(0, neq):
    
    # set maximum distance window
    max_dist = 10**((catalogue.data['magnitude'][i]-4.317)*0.6) \
                     + 17.0/np.sqrt(catalogue.data['magnitude'][i])
    
    # set time-dependent distance cut-off                 
    if (evdate[i] <= test_day_1 ):
        max_dist = max_dist + 5.0
    elif (evdate[i] <= test_day_2 ):
        max_dist = max_dist + 10.0
           
    # for subsequent earthquakes, check distance from last event
    inter_evdist = haversine(catalogue.data['longitude'][i+1:],
                       catalogue.data['latitude'][i+1:],
                       catalogue.data['longitude'][i],
                       catalogue.data['latitude'][i])
                       
    inter_evdist = inter_evdist.flatten()
    
    # get inter-event time
    inter_evtime = evdate[i+1:] - evdate[i]
    inter_evdays = []
    for t in inter_evtime:
        inter_evdays.append(t.days)
    
    # get interevent magnitude
    inter_evmag = delta_mag*catalogue.data['magnitude'][i] - catalogue.data['magnitude'][i+1:]
                       
    #  Flagging aftershocks so test forward in time
    #event_window = evdate[i] + dt.timedelta(days=max_time[i])  
        
    idx = np.where((inter_evdist < max_dist) & (inter_evdays < max_time[i]) \
                    & (inter_evmag > 0.0))[0]
                   
    # set aftershocks
    flagvector[i+1+idx] = 1

    
  
      
      
      
      
      
      
      
      
      
      