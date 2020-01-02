# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:42:50 2017

@author: u56903
"""
def checkstr(num):
    '''
    check if number is nan.  If True, return blank ('')
    '''
    from numpy import isnan
    
    if  isinstance(num, basestring):
        return num
    else:
        if isnan(num):
            return ''
        else:
            return str(num)


def ggcat2ascii(ggcat_dict, outfile):
    '''
    # takes event dictionary format as parsed by catalogues.parsers.parse_ggcat
    
    # exports to outfile
    '''
    # make file for writing
    cattxt = ''
    
    for ev in ggcat_dict:
        # set string constant width
        datestr = '{0.year:4d} {0.month:02d} {0.day:02d} {0.hour:02d}{0.minute:02d}'.format(ev['datetime'])
        
        # get dep format
        depStr = str("%0.1f" % ev['dep'])
        if len(depStr) == 3:
            depStr = '  ' + depStr
        elif len(depStr) == 4:
            depStr = ' ' + depStr
        
        # make line
        line = ' '.join((datestr, str("%0.3f" % ev['lon']), str("%0.3f" % ev['lat']), \
                         depStr, str("%0.2f" % ev['prefmag']), \
                         ev['prefmagtype'], ev['auth']))
                       
        cattxt = cattxt + line + '\n'
                
    #f = open(path.join(outfolder, outfile), 'wb')
    f = open(outfile, 'w')
    f.write(cattxt)
    f.close()
    
def ggcat2hmtk_csv(ggcat_dict, hmtkfile, prefmag):
    '''
    prefmag = 'orig' or 'mw'
    '''
    
    '''
    takes catalogue dictionary format as parsed by catalogues.parsers.parse_ggcat
    
    returns OQ compliant catalogue in csv fmt
    '''
    
    # make oq cat dict
    try:
        test = ggcat_dict[0]['ml2mw_qd']
        header = ','.join(('eventID','year', 'month', 'day', 'hour', 'minute', 'second', \
                           'longitude', 'latitude','depth','magnitude','magnitudeType', \
                           'Agency', 'flag', 'mx_origML', 'mx_origType', 'mx_revML', 'pref_mw', \
                           'ml2mw_qd','ml2mw_bl'))
    except:
        header = ','.join(('eventID','year', 'month', 'day', 'hour', 'minute', 'second', \
                           'longitude', 'latitude','depth','magnitude','magnitudeType', \
                           'Agency', 'flag', 'mx_origML', 'mx_origType', 'mx_revML', 'pref_mw'))
                           
    oq_dat = header + '\n'
    
    '''
    tmpdict = {'auth':line[7], 'place':line[30],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                   'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                   'prefmag':float(line[28]), 'prefmagtype':line[29], 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                   'mw':float(line[8]), 'mp':float(line[17]), 'fixdep':0, 'datetime':evdt, 'dependence':str(line[3]), 'mx_orig':float(line[20]), \
                   'mx_origType':omt, 'mx_rev_ml':float(line[21]), 'mx_rev_src':line[22], 'mw_src':line[-2]}
    '''
                   
                   
    # loop thru eqs
    for ggc in ggcat_dict:
        # make datstr - strftime does not work for dats < 1900!
        datestr = '{0.year:4d}{0.month:02d}{0.day:02d}{0.hour:02d}{0.minute:02d}'.format(ggc['datetime'])
        
        # flag dependent or man-made events
        '''
        # for 2012 catalogue
        if ggc['dependence'] == 'Aftershock' or ggc['dependence'] == 'Foreshock':
            flag = '1'
        elif  ggc['ev_type'] == 'blast' or ggc['ev_type'] == 'coal':
            flag = '2'
        else:
            flag = '0'
        '''
        # for 2018 catalogue
        try:
           if ggc['dependence'] == 0:
               flag = '1'
           else:
               flag = '0'
        except:
            flag = '-99'
            
        # set Original magnitude as main magnitude for declustering (mx_orig) and add additional columns
        if prefmag == 'orig':
            # try alt ml2mw version
            try:
                line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                             checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['mx_origML']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_origML']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_revML']), checkstr(ggc['prefmag']), \
                             checkstr(ggc['ml2mw_qd']),checkstr(ggc['ml2mw_bl'])))
                             
            except:
                line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                             checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['mx_origML']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_origML']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_revML']), checkstr(ggc['prefmag'])))
        
        # else use pref mw
        else:
            try:
                line = ','.join((datestr, checkstr(ggc['datetime'].year), checkstr(ggc['datetime'].month),checkstr(ggc['datetime'].day), \
                             checkstr(ggc['datetime'].hour),checkstr(ggc['datetime'].minute),checkstr(ggc['datetime'].second),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['prefmag']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_origML']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_revML']), checkstr(ggc['prefmag']), \
                             checkstr(ggc['ml2mw_qd']),checkstr(ggc['ml2mw_bl'])))
                                             
            except:
                line = ','.join((datestr, checkstr(ggc['datetime'].year), checkstr(ggc['datetime'].month),checkstr(ggc['datetime'].day), \
                             checkstr(ggc['datetime'].hour),checkstr(ggc['datetime'].minute),checkstr(ggc['datetime'].second),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['prefmag']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_origML']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_revML']), checkstr(ggc['prefmag'])))     
        '''
        # for making MX_REV_ML file        
        line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                         checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                         checkstr(ggc['dep']),checkstr(ggc['mx_rev_ml']),checkstr(ggc['mx_rev_src']),ggc['auth'], flag))
        '''                 
        oq_dat += line + '\n'
        
    #write to OQ out
    print('Writing HMTK csv...')
    f = open(hmtkfile, 'w')
    f.write(oq_dat)
    f.close()
    
def iscgem2hmtk_csv(iscgem_dict, hmtkfile):
    
    '''
    takes catalogue dictionary format as parsed by catalogues.parsers.parse_iscgem
    
    returns OQ compliant catalogue in csv fmt
    '''
    
    # make oq cat dict
    header = ','.join(('eventID','year', 'month', 'day', 'hour', 'minute', 'second', \
                       'longitude', 'latitude','depth','magnitude','magnitudeType', \
                       'Agency', 'flag'))

    oq_dat = header + '\n'
                       
    for ev in iscgem_dict:
        line = ','.join((ev['eventid'], str(ev['datetime'].year), str(ev['datetime'].month), str(ev['datetime'].day), \
                         str(ev['datetime'].hour), str(ev['datetime'].minute), str(ev['datetime'].second), str(ev['lon']), \
                         str(ev['lat']), str(ev['dep']), str(ev['mw']), 'MW', ev['mo_auth'], '0'))
                         
        oq_dat += line + '\n'
        
    #write to OQ out
    print('Writing HMTK csv...')
    f = open(hmtkfile, 'w')
    f.write(oq_dat)
    f.close()   

# writes htmk dict to shapefile
def htmk2shp(cat, outshp):
    '''
    cat = htmk dictionory
    outshp = output shapefile
    '''
    import shapefile
    from numpy import isnan
    
    print('Making shapefile...')
    w = shapefile.Writer(shapefile.POINT)
    w.field('EVID','C','15')
    w.field('AGENCY','C','15')
    w.field('YEAR','F', 4, 0)
    w.field('MONTH','F', 2, 0)
    w.field('DAY','F', 2, 0)
    w.field('HOUR','F', 2, 0)
    w.field('MIN','F', 2, 0)
    w.field('SEC','F', 2, 0)
    w.field('LON','F', 10, 4)
    w.field('LAT','F', 10, 4)    
    w.field('DEP','F', 10, 2)
    w.field('PREF_MW','F', 5, 2)
    w.field('ORIG_MAG','F', 5, 2)
    w.field('ORIG_MAG_TYPE','C','6')
    
    # now loop thru records
    for i in range(0, len(cat.data['eventID'])):
        if isnan(cat.data['longitude'][i]) | isnan(cat.data['latitude'][i]):
            lon = 0.0
            lat = 0.0
        else:
            lon = round(cat.data['longitude'][i],4)
            lat = round(cat.data['latitude'][i],4)
    
        w.point(lon, lat)
        w.record(cat.data['eventID'][i],cat.data['Agency'][i], \
                 cat.data['year'][i],cat.data['month'][i], \
                 cat.data['day'][i],cat.data['hour'][i], \
                 cat.data['minute'][i],cat.data['second'][i], \
                 round(cat.data['longitude'][i],4), \
                 round(cat.data['latitude'][i],4), \
                 cat.data['depth'][i],round(cat.data['pref_mw'][i],2),\
                 cat.data['mx_origML'][i],cat.data['mx_origType'][i])
    
    print('Writing shapefile...')
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'w')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()
    
# writes htmk dict to shapefile
def htmk2shp_isc(cat, outshp):
    '''
    cat = htmk dictionory
    outshp = output shapefile
    '''
    import shapefile
    from numpy import isnan
    
    print('Making shapefile...')
    w = shapefile.Writer(shapefile.POINT)
    w.field('EVID','C','15')
    w.field('AGENCY','C','15')
    w.field('YEAR','F', 4, 0)
    w.field('MONTH','F', 2, 0)
    w.field('DAY','F', 2, 0)
    w.field('HOUR','F', 2, 0)
    w.field('MIN','F', 2, 0)
    w.field('SEC','F', 2, 0)
    w.field('LON','F', 10, 4)
    w.field('LAT','F', 10, 4)    
    w.field('DEP','F', 10, 2)
    w.field('MAGNITUDE','F', 6, 2)
    
    # now loop thru records
    for i in range(0, len(cat.data['eventID'])):
        if isnan(cat.data['longitude'][i]) | isnan(cat.data['latitude'][i]):
            lon = 0.0
            lat = 0.0
        else:
            lon = round(cat.data['longitude'][i],4)
            lat = round(cat.data['latitude'][i],4)
    
        w.point(lon, lat)
        w.record(cat.data['eventID'][i],cat.data['Agency'][i], \
                 cat.data['year'][i],cat.data['month'][i], \
                 cat.data['day'][i],cat.data['hour'][i], \
                 cat.data['minute'][i],cat.data['second'][i], \
                 round(cat.data['longitude'][i],4), \
                 round(cat.data['latitude'][i],4), \
                 cat.data['depth'][i], cat.data['magnitude'][i])
    
    print('Writing shapefile...')
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'w')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()

    
    
    
    
    
    
    
    
    
    
    
    
    