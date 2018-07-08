
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
        
        # get datetime
        try:
            evdt = datetime.strptime(line[0], '%Y-%m-%d %H:%M:%S')
        except:
            evdt = datetime.strptime(line[0], '%Y-%m-%d %H:%M')
        
        # get original magnitude type
        omt = str(line[22]).strip('REV_')
        
            
        # fill temp dict
        try:        
            # V0.12
            tmpdict = {'auth':line[7], 'place':line[30],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                       'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                       'prefmag':float(line[28]), 'prefmagtype':line[29], 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                       'mw':float(line[8]), 'mp':float(line[17]), 'fixdep':0, 'datetime':evdt, 'dependence':str(line[3]), 'mx_orig':float(line[20]), \
                       'mx_origType':omt, 'mx_rev_ml':float(line[21]), 'mx_rev_src':line[22], 'mw_src':line[-2], 'ev_type':str(line[2])}
        except:
            # V0.11
            tmpdict = {'auth':line[7], 'place':line[29],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                       'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                       'prefmagtype':line[28], 'prefmag':float(line[27]), 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                       'mw':float(line[8]), 'mp':float(line[17]), 'fixdep':0, 'datetime':evdt, 'dependence':str(line[3]), 'mx_orig':float(line[20]), \
                       'mx_rev_ml':float(line[21]), 'mx_rev_src':line[22]}
        
        austcat.append(tmpdict)
        
    return austcat
    
# parses catalogue of format used for the NSHA2018
def parse_NSHA2018_catalogue(nsha2018cat):
    '''
    nsha2018cat: catalogue is csv format
    '''
    import csv
    from numpy import nan
    from datetime import datetime
    #from os import path
    
    # for testing only
    #nsha2018cat = path.join('data', 'NSHA18CAT.MW.V0.1.csv')
    
    # open file
    raw = open(nsha2018cat).readlines()[1:] # exclude header
    
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
        
        # get datetime
        try:
            evdt = datetime.strptime(line[0], '%Y-%m-%d %H:%M:%S')
        except:
            evdt = datetime.strptime(line[0], '%Y-%m-%d %H:%M')
        
        # read V0.2 with alternatle ML mag conversion
        try:
            '''
            ml2mw_qds = quadratic simulated
            ml2mw_qde = quadratic empirical
            ml2mw_ble = bi-linear empirical
            '''
            # fill temp dict
            tmpdict = {'auth':line[7], 'place':line[29],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                       'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                       'prefmag':float(line[28]), 'prefmagtype':line[29], 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                       'mw':float(line[8]), 'fixdep':0, 'datetime':evdt, 'dependence':int(line[3]), 'mx_origML':float(line[18]), \
                       'mx_origType':str(line[19]), 'mx_revML':float(line[20]), 'mx_rev_src':line[22], 'mw_src':line[28], 'ev_type':str(line[2]), \
                       'ml2mw_qds':float(line[25]),'ml2mw_ble':float(line[26]),'ml2mw_qde':float(line[27])}
        
        # read V0.1
        except:
            # fill temp dict
            tmpdict = {'auth':line[7].strip(), 'place':line[28],'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                       'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[4]), 'lat':float(line[5]), 'dep':float(line[6]), \
                       'prefmag':float(line[26]), 'prefmagtype':line[27], 'ml':float(line[14]), 'mb':float(line[12]), 'ms':float(line[10]), \
                       'mw':float(line[8]), 'fixdep':0, 'datetime':evdt, 'dependence':int(line[3]), 'mx_origML':float(line[18]), \
                       'mx_origType':str(line[19]), 'mx_revML':float(line[20]), 'mx_rev_src':line[22], 'mw_src':line[27], 'ev_type':str(line[2]), \
                       'ml_region':str(line[16])}
    
        austcat.append(tmpdict)
            
    return austcat
    
def parse_altmag_hmtk_catalogue(hmtk_csv):
    from datetime import datetime
    
    lines = open(hmtk_csv).readlines()[1:]
    
    # set arrays
    altMWdict = []
    
    # fill arrays
    for line in lines:
        data = line.strip().split(',')
        dateStr = data[0]
        lon = float(data[9])
        lat = float(data[10])
        dep = float(data[14])
        mx_origML = float(data[19])
        mx_revML = float(data[21])
        mw_pref = float(data[22])
        mw_qds = float(data[25])
        mw_ble = float(data[23])
        mw_qde = float(data[24])
        
        # set empirical alt MW
        if mw_qds == mw_pref:
            mw_alt_ble = mw_ble
            mw_alt_qde = mw_qde
        else:
            mw_alt_ble = mw_pref
            mw_alt_qde = mw_pref
            
        # make datetime object
        try:
            evdt = datetime.strptime(dateStr, '%Y%m%d%H%M')
        except:
            print dateStr
            evdt = datetime.strptime(dateStr, '%Y-%m-%d %H:%M')
        ev_date = evdt
        
        tmpdict = {'datetime':ev_date, 'lon':lon, 'lat':lat, 'dep':dep,
                   'mx_origML':mx_origML, 'mx_revML': mx_revML,
                   'mw_pref':mw_pref, 'mw_qds':mw_qds, 'mw_qde':mw_qde, 
                   'mw_ble':mw_ble, 'mw_alt_ble':mw_alt_ble, 
                   'mw_alt_qde':mw_alt_qde, 'prefmag':mw_pref}
                   	
        altMWdict.append(tmpdict)
        
    return altMWdict, len(altMWdict)


def parse_iscgem(iscgemcsv):
    """
    function to parse the ISC-GEM V5 earthquake catalogue in csv format
    
    returns a list of dictionaries, each dictionary relating to a single event
    """
    
    import csv
    from numpy import nan, isnan, floor
    from catalogue.parsers import checkint
    from datetime import datetime
    
    # open file
    raw = open(iscgemcsv).readlines()[60:] # exclude header    
    
    # parse csv    
    lines = csv.reader(raw)
    
    # set array to append event dictionaries
    iscgemCat = []
    
    for line in lines:
        evdt = datetime.strptime(line[0].strip(), '%Y-%m-%d %H:%M:%S.%f')
        tmpdict = {'eventid':line[-1].strip(), 'datetime':evdt, 'year':evdt.year, 'month':evdt.month, 'day':evdt.day, \
                   'hour':evdt.hour, 'min':evdt.minute, 'sec':evdt.second, 'lon':float(line[2].strip()), 'lat':float(line[1].strip()), \
                   'dep':float(line[7].strip()), 'mw':float(line[10].strip()), 'mo_auth':line[16].strip()}
                   
        iscgemCat.append(tmpdict)
            
    return iscgemCat
    
    
    
    
    
    
    
    
    
    
    