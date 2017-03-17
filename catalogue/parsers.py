
# checks if string can be converted to an int
# returns nan if not
def checkint(intstr):
    try:
        return int(intstr)
    except:
        from numpy import nan
        return nan

def parse_ggcat(ggcatcsv):
    """
    function to parse Gary Gibson's earthquake catalogue in csv format
    
    returns a list of dictionaries, each dictionary relating to a single event
    """
    
    import csv
    from numpy import nan, isnan, floor
    from catalogue.parsers import checkint
    import datetime as dt
    
    # open file
    raw = open(ggcatcsv).readlines()[1:] # exclude header    
    
    # parse csv    
    lines = csv.reader(raw)
    
    # set array to append event dictionaries
    ggcat = []
    
    # loop through events
    for line in lines:
        # set null vals to nans
        for i in range(len(line)):
           if len(str(line[i]))<1:
              line[i] = nan
        
        # fill temp dict
        tmpdict = {'auth':line[0], 'place':line[1],'year':checkint(line[5]), 'month':checkint(line[6]), 'day':checkint(line[7]), \
                   'hour':checkint(line[8]), 'min':checkint(line[9]), 'sec':float(line[10]), 'lon':float(line[11]), 'lat':float(line[12]), 'dep':float(line[13]), \
                   'zcode':line[14], 'prefmagtype':line[15], 'prefmag':float(line[16]), 'ml':float(line[17]), 'mb':float(line[18]), 'ms':float(line[19]), \
                   'mw':float(line[20]), 'md':float(line[21]), 'mp':float(line[22]), 'fixdep':0}
                   	
        # add datetime
        if ~isnan(tmpdict['sec']):
            if int(floor(tmpdict['sec'])) >= 60:
                tmpdict['sec'] = 59
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'], tmpdict['hour'], tmpdict['min'], int(floor(tmpdict['sec'])))
        elif ~isnan(tmpdict['hour']):
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'], tmpdict['hour'], tmpdict['min'], 0)
        else:
            tmpdict['datetime'] = dt.datetime(tmpdict['year'], tmpdict['month'], tmpdict['day'])
        
        ggcat.append(tmpdict)
        
    return ggcat

# parses catalogue of format used for the NSHA2012
def parse_NSHA2012_catalogue(nsha2012cat):
    '''
    nsha2012cat: catalogue is csv format
    '''
    import csv
    from numpy import nan
    from datetime import datetime
    
    # for testing only
    #nsha2012cat = path.join('data', 'AUSTCAT.MW.V0.11.csv')
    
    # open file
    raw = open(nsha2012cat).readlines()[1:] # exclude header    
    
    # parse csv    
    lines = csv.reader(raw)
    
    # set array to append event dictionaries
    austcat = []
    
    # loop through events
    for line in lines:
        # set null vals to nans
        for i in range(len(line)):
           if len(str(line[i]))<1:
              line[i] = nan
        '''
        line = ','.join((datestr, checkstr(ggc['year']), checkstr(ggc['month']),checkstr(ggc['day']), \
                             checkstr(ggc['hour']).zfill(2),checkstr(ggc['min']).zfill(2),checkstr(ggc['sec']),checkstr(ggc['lon']),checkstr(ggc['lat']), \
                             checkstr(ggc['dep']),checkstr(ggc['prefmag']),ggc['prefmagtype'],ggc['auth']))
        '''
        # get datetime
        evdt = datetime.strptime(line[0], '%Y-%m-%d %H:%M:%S')
            
        # fill temp dict
        tmpdict = {'auth':line[7], 'place':line[29],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                   'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                   'prefmagtype':line[28], 'prefmag':float(line[27]), 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                   'mw':float(line[8]), 'mp':float(line[17]), 'fixdep':0, 'datetime':evdt, 'dependence':str(line[3])}
        
        austcat.append(tmpdict)
        
    return austcat
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        