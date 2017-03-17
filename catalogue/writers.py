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
    
def ggcat2hmtk_csv(ggcat_dict, hmtkfile):
    
    '''
    takes catalogue dictionary format as parsed by catalogues.parsers.parse_ggcat
    
    returns OQ compliant catalogue in csv fmt
    '''
    
    # make oq cat dict
    header = ','.join(('eventID','year', 'month', 'day', 'hour', 'minute', 'second', 'longitude', 'latitude','depth','magnitude','magnitudeType','Agency', 'flag'))
    oq_dat = header + '\n'
    
    # loop thru eqs
    for ggc in ggcat_dict:
        # make datstr - strftime does not work for dats < 1900!
        datestr = '{0.year:4d}{0.month:02d}{0.day:02d}{0.hour:02d}{0.minute:02d}'.format(ggc['datetime'])
        
        if ggc['dependence'] == 'Aftershock' or ggc['dependence'] == 'Foreshock':
            flag = '1'
        else:
            flag = '0'
        
        line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                         checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                         checkstr(ggc['dep']),checkstr(ggc['prefmag']),checkstr(ggc['prefmagtype']),ggc['auth'], flag))
        oq_dat += line + '\n'
        
    #write to OQ out
    print 'Writing HMTK csv...'
    f = open(hmtkfile, 'wb')
    f.write(oq_dat)
    f.close()