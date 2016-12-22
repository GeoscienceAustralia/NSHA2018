

def parse_ggcat(ggcatcsv):
    """
    function to parse Gary Gibson's earthquake catalogue in csv format
    
    returns a list of dictionaries, each dictionary relating to a single event
    """
    
    import csv
    from numpy import nan, isnan, floor
    from misc_tools import checkint
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