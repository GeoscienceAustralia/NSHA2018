# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 13:30:27 2017

@author: u56903
"""

from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from hmtk.seismicity.utils import haversine
from parsers import parse_NSHA2018_catalogue
from writers import ggcat2hmtk_csv, htmk2shp
import numpy as np
import datetime as dt
from os import path, remove
from copy import deepcopy

# does the grunt work to decluster catalogue
def flag_dependent_events(catalogue, flagvector, doAftershocks, method):
	
    '''
    catalogue: dictionary of earthquakes in HMTK catalogue format, 
               parsed using CsvCatalogueParser
    flagvector: integer vector of length of catalogue
    doAftershocks: 
        if == True: decluster aftershocks
        if == False: decluster foreshocks
    method: either "Leonard08" or "Stein08"
    '''
    
    # get number of events
    neq = len(catalogue.data['magnitude'])
    
    # set periods of confidence
    test_day_1 = dt.datetime(1960,1,1) # Earthquakes older than this are assumed to be very poorly located.
    test_day_2 = dt.datetime(1970,1,1) # Earthquakes older than this are assumed to be poorly located.
    
    # set delta-magnitude threshold
    delta_mag = 0.989 # Aftershocks must be less than 0.989 m.u. of the main shock (95% of the moment).
    
    # set time window
    if method == 'Leonard08':
        max_time_foreshock = 10**((catalogue.data['magnitude']-1.85)*0.69)
        max_time_aftershock = 10**((catalogue.data['magnitude']-2.70)*1.1) + 4.0
    elif method == 'Stien08':
        max_time_foreshock = 10**((catalogue.data['magnitude']-2.70)*1.1) + 4.0
        max_time_aftershock = 10**((catalogue.data['magnitude']-2.70)*1.1) + 4.0
        
    # get event time datevector
    evdate = []
    for i in range(0, neq):
    
        # get event datetime
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
        # flag aftershocks
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
            inter_evmag = delta_mag*catalogue.data['magnitude'][i] \
                          - catalogue.data['magnitude'][i+1:]
                               
            # now find aftershocks to flag
            idx = np.where((inter_evdist < max_dist) & (inter_evdays < max_time_aftershock[i]) \
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
            inter_evmag = delta_mag*catalogue.data['magnitude'][i] \
                          - catalogue.data['magnitude'][0:i]
                               
            # now find aftershocks to flag
            idx = np.where((inter_evdist < max_dist) & (inter_evdays < max_time_foreshock[i]) \
                            & (inter_evmag > 0.0))[0]
                            
            # set foreshock flag
            flagvector[idx] = 1
            
    
    return flagvector

# methods to call Leonard 2008 & Stein 2008 declustering algorthms 
def decluster_SCR(method, cat, deblastOnly):

    # set flag for dependent events
    flagvector = np.zeros(len(cat.data['magnitude']), dtype=int)
    
    #########################################################################
    # call declustering
    #########################################################################
    # use "method" to decluster
    if deblastOnly == False:
        # flag aftershocks
        doAftershocks = True
        flagvector_as = flag_dependent_events(cat, flagvector, doAftershocks, method)
        
        # flag foreshocks
        doAftershocks = False
        flagvector_asfs = flag_dependent_events(cat, flagvector_as, doAftershocks, method)
        
        # now find mannually picked foreshocks/aftershocks (1 = fore/aftershock; 2 = blast/coal)
        # idx = np.where(cat.data['flag'] >= 1)[0] # for 2012 version
        
        # now find mannually picked foreshocks/aftershocks (1 = dependent events)
        # !!! REMEMBER - "writers.ggcat2hmtk_csv" sets dependent events flag to 1 !!!
        idx = np.where(cat.data['flag'] == 1)[0]
        flagvector_asfsman = flagvector_asfs
        flagvector_asfsman[idx] = 1
    
    # else remove coal & blast events only
    else:
        idx = np.where(cat.data['flag'] == 0)[0]
        flagvector_asfsman = flagvector
        flagvector_asfsman[idx] = 1
            
    
    #########################################################################
    # purge non-poissonian events
    #########################################################################
    
    # adding cluster flag to the catalog
    cat.data['cluster_flag'] = flagvector_asfsman
    
    # create a copy from the catalogue object to preserve it
    catalogue_l08 = deepcopy(cat)
    
    catalogue_l08.purge_catalogue(flagvector_asfsman == 0) # cluster_flags == 0: mainshocks
    
    print 'Leonard 2008\tbefore: ', cat.get_number_events(), " after: ", catalogue_l08.get_number_events()
    
    #####################################################
    # write declustered catalogue
    #####################################################
    
    # setup the writer
    declustered_catalog_file = '.'.join(hmtk_csv.split('.')[0:2])+'_declustered_test.csv'
    
    # if it exists, delete previous file
    try:
        remove(declustered_catalog_file)
    except:
        print declustered_catalog_file, 'does not exist' 
    
    # set-up writer
    writer = CsvCatalogueWriter(declustered_catalog_file) 
    
    # write
    writer.write_file(catalogue_l08)
    
    print 'Declustered catalogue: ok\n'
    
    return declustered_catalog_file, catalogue_l08


# do Gardner & Knopoff (1974 declustering)
def decluster_GK74(catalogue):
    from hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1
    from hmtk.seismicity.declusterer.distance_time_windows import GardnerKnopoffWindow
    
    decluster_config = {'time_distance_window': GardnerKnopoffWindow(),
                        'fs_time_prop': 1.0}
    
    #####################################################
    # decluster here
    #####################################################
    
    print 'Running GK declustering...'
    decluster_method = GardnerKnopoffType1()
    
    #---------------------------------------------
    # declustering
    cluster_index_gk, cluster_flags_gk = decluster_method.decluster(catalogue, decluster_config)
    #---------------------------------------------
    
    # adding to the catalog
    # The cluster flag (main shock or after/foreshock) and cluster index to the catalogue keys
    catalogue.data['cluster_index_gk'] = cluster_index_gk
    catalogue.data['cluster_flags_gk'] = cluster_flags_gk
    
    #####################################################
    # purge remove non-poissonian events
    #####################################################
    
    # create a copy from the catalogue object to preserve it
    catalogue_gk = deepcopy(catalogue)
    
    catalogue_gk.purge_catalogue(cluster_flags_gk == 0) # cluster_flags == 0: mainshocks
    
    print 'Gardner-Knopoff\tbefore: ', catalogue.get_number_events(), " after: ", catalogue_gk.get_number_events()
    
    #####################################################
    # write declustered catalogue
    #####################################################
    
    # setup the writer
    declustered_catalog_file = hmtk_csv.split('.')[0]+'_GK74_declustered.csv'
    
    # if it exists, delete previous file
    try:
        remove(declustered_catalog_file)
    except:
        print declustered_catalog_file, 'does not exist'
    
    # set-up writer
    writer = CsvCatalogueWriter(declustered_catalog_file) 
    
    # write
    writer.write_file(catalogue_gk)
    #writer.write_file(catalogue_af)
    print 'Declustered catalogue: ok\n'


#########################################################################
# !!!!start main code here!!!!
#########################################################################

leonard = True
deblastOnly = False # remove blasts and coal events - if True, does not decluster

#########################################################################
# parse calalogue & convert to HMTK
#########################################################################
'''
recommend declustering on original magnitudes given the declustering algorithm 
was based on these magnitudes
'''
prefmag1 = 'orig' # declusters based on original catalogue magnitude

# Use 2018 NSHA catalogue
print 'Declustering NSHA18CAT.MW.V0.1.csv'
nsha2018csv = path.join('data', 'NSHA18CAT.MW.V0.2.csv')
nsha_dict = parse_NSHA2018_catalogue(nsha2018csv)

# set HMTK file name
if prefmag1 == 'orig':
    hmtk_csv = nsha2018csv.split('.')[0] + '_V0.1_hmtk_mx_orig.csv'
    #hmtk_csv = nsha2018csv.split('.')[0] + '_V0.12_hmtk_mp_orig.csv'
elif prefmag1 == 'mw':
    hmtk_csv = nsha2018csv.split('.')[0] + '_V0.1_hmtk.csv'

# write HMTK csv
ggcat2hmtk_csv(nsha_dict, hmtk_csv, prefmag1)

# parse HMTK csv
parser = CsvCatalogueParser(hmtk_csv)
nshacat = parser.read_file()

cat = nshacat

# make shapefile of full catalogue
htmk2shp(cat, path.join('shapefiles', 'NSHA18CAT_full.shp')) # for testing

#########################################################################
# if Leonard == True
#########################################################################

if leonard == True:
    # try following HMTK format
    method = 'Leonard08'
    if method == 'Leonard08':
        print  'Using Leonard 2008 method...'
    if method == 'Stein08':
        print  'Using Stein 2008 method...'
        
    declustered_catalog_filename, declustered_cat = decluster_SCR(method, cat, deblastOnly)
    
    # make shapefile of declustered catalogue
    htmk2shp(declustered_cat, path.join('shapefiles', 'NSHA18CAT_declustered.shp'))

# do GK74    
else: 
    decluster_GK74(cat)
   
#########################################################################
# merge extra columns removed by HMTK parser
#########################################################################
print 'Merging stripped columns...\n'

prefmag2 = 'mw' # replaces orig mag with preferred MW in declustered catalogue
#prefmag2 = 'orig' # keeps orig mag in declustered catalogue

from misc_tools import checkfloat, dictlist2array

# for testing
#declustered_catalog_filename = 'data/AUSTCAT_V0.12_hmtk_mx_orig_declustered_test.csv'

# get data arrays from original catalogue
for key in nsha_dict[0].keys():
    exec(key + ' = dictlist2array(nsha_dict, key)')

# make datestr list
datestr = []
for i in range(0, len(nsha_dict)):
    datestr.append(''.join((str(nsha_dict[i]['year']), str('%02d' % nsha_dict[i]['month']), \
                       str('%02d' % nsha_dict[i]['day']), str('%02d' % nsha_dict[i]['hour']), \
                       str('%02d' % nsha_dict[i]['min']))))

datestr = np.array(datestr)

# parse declustered file and add cols
lines = open(path.join(declustered_catalog_filename)).readlines()

# set new header
newheader = lines[0].strip() #+ ',mx_origML,mx_origType,mx_revML,pref_mw\n'
newtxt = newheader

#  all new variables are invoked using the "exec" command above
for line in lines[1:]:
    dat = line.strip().split(',')    
    
    # match earthquake info
    if np.isnan(float(dat[16])):
        eqidx = np.where((datestr == dat[0]) \
                         & (lon >= checkfloat(dat[9])-0.01) & (lon <= checkfloat(dat[9])+0.01) \
                         & (lat >= checkfloat(dat[10])-0.01) & (lat <= checkfloat(dat[10])+0.01))[0]
    else:
        eqidx = np.where((datestr == dat[0]) \
                         & (lon >= checkfloat(dat[9])-0.01) & (lon <= checkfloat(dat[9])+0.01) \
                         & (lat >= checkfloat(dat[10])-0.01) & (lat <= checkfloat(dat[10])+0.01) \
                         & (mx_orig >= checkfloat(dat[16])-0.1) & (mx_orig <= checkfloat(dat[16])+0.1))[0]
    
    # testing
    if len(eqidx) > 1:
        print 'Possible duplicate:', datestr[eqidx[0]], eqidx
    
    # replace orig mag in "magnitude" column
    if len(eqidx) > 0:
        # prefmag = MW
        if prefmag2 == 'mw':
            newline = ','.join(dat[0:16]) + ',' + str('%0.2f' % prefmag[eqidx][0]) + ',,MW,' \
                      + ','.join((str('%0.2f' % mx_orig[eqidx][0]), mx_origType[eqidx][0], \
                                  str('%0.2f' % mx_rev_ml[eqidx][0]), str('%0.2f' % prefmag[eqidx][0]))) + '\n'
                                
        # keep orig mag in "magnitude" column
        else:
            newline = ','.join((line.strip(), str('%0.2f' % mx_orig[eqidx][0]), mx_origType[eqidx][0], \
                                str('%0.2f' % mx_rev_ml[eqidx][0]), str('%0.2f' % prefmag[eqidx][0]))) + '\n'
        
        newtxt += newline
    else:
        print 'No match:', line

# re-write declustered catalogue
if deblastOnly == False:
    if prefmag2 == 'mw':
        declustered_csv = nsha2018csv.split('.')[0] + '_V0.1_hmtk_declustered_test.csv'   
    else:
        declustered_csv = declustered_catalog_filename
else:
    if prefmag2 == 'mw':
        declustered_csv = nsha2018csv.split('.')[0] + '_V0.12_hmtk_deblast.csv'

f = open(declustered_csv, 'wb')
f.write(newtxt)
f.close()
  

