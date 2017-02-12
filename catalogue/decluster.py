# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 13:30:27 2017

@author: u56903
"""

from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from hmtk.seismicity.utils import decimal_year, haversine
from parsers import parse_ggcat
from writers import ggcat2hmtk_csv
import numpy as np
import datetime as dt
from os import path, remove
from copy import deepcopy

# does the grunt work to decluster catalogue
def flag_dependent_events(catalogue, flagvector, doAftershocks, method):
	
    '''
    catalogue: dictionary of earthquakes in HMTK catalogue format
    flagvector: integer vector of length of catalogue
    doAftershocks: 
        if == True: decluster aftershocks
        if == False: decluster foreshocks
    method: either "Leonard08" of "Stein08"
    '''
    
    # get number of events
    neq = len(catalogue.data['magnitude'])  # Number of earthquakes
    
    # set periods of confidence
    test_day_1 = dt.datetime(1960,1,1) # Earthquakes older than this are assumed to be very poorly located.
    test_day_2 = dt.datetime(1970,1,1) # Earthquakes older than this are assumed to be poorly located.
    
    # set delta-magnitude threshold
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

    if doAftershocks == True:
        print 'Flagging aftershocks...'
    else:
        print 'Flagging foreshocks...'
    
    # loop through earthquakes
    for i in range(0, neq):
        
        # set maximum distance window
        max_dist = 10**((catalogue.data['magnitude'][i]-4.317)*0.6) \
                         + 17.0/np.sqrt(catalogue.data['magnitude'][i])
        
        # set time-dependent distance cut-off
        if (evdate[i] <= test_day_1 ):
            max_dist = max_dist + 5.0
        elif (evdate[i] <= test_day_2 ):
            max_dist = max_dist + 10.0
        
        #########################################################################
        # flag foreshocks
        #########################################################################
        
        if doAftershocks == True:
        
            # for subsequent earthquakes, check distance from current event
            inter_evdist = haversine(catalogue.data['longitude'][i+1:],
                                     catalogue.data['latitude'][i+1:],
                                     catalogue.data['longitude'][i],
                                     catalogue.data['latitude'][i])
             
             # flatten distance array
            inter_evdist = inter_evdist.flatten()
            
            # get inter-event time in days
            inter_evtime = evdate[i+1:] - evdate[i]
            inter_evdays = []
            for t in inter_evtime:
                inter_evdays.append(t.days)
            
            # get interevent magnitude
            inter_evmag = delta_mag*catalogue.data['magnitude'][i] - catalogue.data['magnitude'][i+1:]
                               
            # now find aftershocks to flag
            idx = np.where((inter_evdist < max_dist) & (inter_evdays < max_time[i]) \
                            & (inter_evmag > 0.0))[0]
                            
            # set aftershock flag
            flagvector[i+1+idx] = 1
            
        #########################################################################
        # flag foreshocks
        #########################################################################

        elif doAftershocks == False:
        
            # for earlier earthquakes, check distance from current event
            inter_evdist = haversine(catalogue.data['longitude'][0:i],
                                     catalogue.data['latitude'][0:i],
                                     catalogue.data['longitude'][i],
                                     catalogue.data['latitude'][i])
             
             # flatten distance array
            inter_evdist = inter_evdist.flatten()
            
            # get inter-event time in days
            inter_evtime = evdate[i] - evdate[0:i]
            inter_evdays = []
            for t in inter_evtime:
                inter_evdays.append(t.days)
            
            # get interevent magnitude
            inter_evmag = delta_mag*catalogue.data['magnitude'][i] - catalogue.data['magnitude'][0:i]
                               
            # now find aftershocks to flag
            idx = np.where((inter_evdist < max_dist) & (inter_evdays < max_time[i]) \
                            & (inter_evmag > 0.0))[0]
                            
            # set foreshock flag
            flagvector[idx] = 1
        
    return flagvector

# !!!!start main code here!!!!

#########################################################################
# parse calalogue & convert to HMTK
#########################################################################

ggcatcsv = path.join('data', 'GGcat-161025.csv')
ggdict = parse_ggcat(ggcatcsv)

# set HMTK file name
hmtk_csv = ggcatcsv.split('.')[0] + '_hmtk.csv'

# write HMTK csv
ggcat2hmtk_csv(ggdict, hmtk_csv)

# parse HMTK csv
parser = CsvCatalogueParser(hmtk_csv)
ggcat = parser.read_file()

#########################################################################
# set variables for declustering
#########################################################################

# try following HMTK format
method = 'Leonard08'
if method == 'Leonard08':
    print  'Using Leonard 2008 method...'
if method == 'Stein08':
    print  'Using Stein 2008 method...'

# rename catalogue
catalogue = ggcat

# get number of events
#neq = len(catalogue.data['magnitude'])  # Number of earthquakes

# set flag for dependent events
flagvector = np.zeros(len(catalogue.data['magnitude']), dtype=int)

#########################################################################
# call declustering
#########################################################################

# flag aftershocks
doAftershocks = True
flagvector_as = flag_dependent_events(catalogue, flagvector, doAftershocks, method)

# flag aftershocks
doAftershocks = False
flagvector_asfs = flag_dependent_events(catalogue, flagvector_as, doAftershocks, method)

#########################################################################
# purge non-poissonian events
#########################################################################

# adding cluster flag to the catalog
catalogue.data['cluster_flag'] = flagvector_asfs

# create a copy from the catalogue object to preserve it
catalogue_l08 = deepcopy(catalogue)

catalogue_l08.purge_catalogue(flagvector_asfs == 0) # cluster_flags == 0: mainshocks

print 'Leonard 2008\tbefore: ', catalogue.get_number_events(), " after: ", catalogue_l08.get_number_events()

#####################################################
# write declustered catalogue
#####################################################

# setup the writer
declustered_catalog_file = hmtk_csv.split('.')[0]+'_declustered.csv'

# if it exists, delete previous file
remove(declustered_catalog_file)

# set-up writer
writer = CsvCatalogueWriter(declustered_catalog_file) 

# write
writer.write_file(catalogue_l08)

print 'Declustered catalogue: ok\n'
