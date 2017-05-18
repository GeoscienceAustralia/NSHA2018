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
        
        # make line        
        line = ' '.join((datestr, str("%0.3f" % ev['lon']), str("%0.3f" % ev['lat']), \
                         str("%0.1f" % ev['dep']), str("%0.1f" % ev['prefmag']), \
                         ev['prefmagtype'], ev['auth']))
                       
        cattxt = cattxt + line + '\n'
                
    #f = open(path.join(outfolder, outfile), 'wb')
    f = open(outfile, 'wb')
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
        if ggc['dependence'] == 'Aftershock' or ggc['dependence'] == 'Foreshock' \
           or ggc['ev_type'] == 'blast' or ggc['ev_type'] == 'coal':
            flag = '1'
        else:
            flag = '0'
        
        # set Original magnitude as main magnitude for declustering (mx_orig) and add additional columns
        if prefmag == 'orig':
            line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                             checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['mx_orig']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_orig']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_rev_ml']), checkstr(ggc['prefmag'])))
        
        # else use pref mw
        else:
            line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                             checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['prefmag']),checkstr(ggc['mx_origType']),ggc['auth'], flag, \
                             checkstr(ggc['mx_orig']), checkstr(ggc['mx_origType']), checkstr(ggc['mx_rev_ml']), checkstr(ggc['prefmag'])))     
        '''
        # for making MX_REV_ML file        
        line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                         checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                         checkstr(ggc['dep']),checkstr(ggc['mx_rev_ml']),checkstr(ggc['mx_rev_src']),ggc['auth'], flag))
        '''                 
        oq_dat += line + '\n'
        
    #write to OQ out
    print 'Writing HMTK csv...'
    f = open(hmtkfile, 'wb')
    f.write(oq_dat)
    f.close()