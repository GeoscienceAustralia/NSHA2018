'''
Removes duplicate events from 

'''
from parsers import parse_NSHA2012_catalogue
from misc_tools import dictlist2array, checkfloat
from numpy import array, isnan, where, hstack, delete
from os import path

# Parse 2012 NSHA catalogue as dict
nsha2012csv = path.join('data', 'AUSTCAT.MW.V0.12.csv')
nsha2012csv = path.join('data', 'AUSTCAT.MP.V0.12.csv') # for use with preferred magnitudes
nsha2012csv = path.join('data', 'AUSTCAT.MW.V0.11.csv') # for use with 2012 catalogue

# Parse 2012 NSHA catalogue as dict
nsha_dict = parse_NSHA2012_catalogue(nsha2012csv)

# get data arrays from original catalogue
print array(nsha_dict[0].keys()), '\n'
for key in nsha_dict[0].keys():
    exec(key + ' = dictlist2array(nsha_dict, key)')

# make datestr array
datestr = []
for i in range(0, len(nsha_dict)):
    datestr.append(''.join((str(nsha_dict[i]['year']), str('%02d' % nsha_dict[i]['month']), \
                       str('%02d' % nsha_dict[i]['day']), str('%02d' % nsha_dict[i]['hour']), \
                       str('%02d' % nsha_dict[i]['min']))))

datestr = array(datestr)

# parse csv file as lines too - we will copy and paste these!
lines = open(nsha2012csv).readlines()

# set new next & header
newtxt = lines[0]

# loop thru records and find duplicates
i = 0
delidx = []
inidx  = []
while i < len(datestr):
    
    # match earthquake info
    if isnan(mx_orig[i]):
        eqidx = where((datestr == datestr[i]) \
                      & (lon >= lon[i]-0.01) & (lon <= lon[i]+0.01) \
                      & (lat >= lat[i]-0.01) & (lat <= lat[i]+0.01))[0]
    else:
        eqidx = where((datestr == datestr[i]) \
                      & (lon >= lon[i]-0.01) & (lon <= lon[i]+0.01) \
                      & (lat >= lat[i]-0.01) & (lat <= lat[i]+0.01) \
                      & (mx_orig >= mx_orig[i]-0.1) & (mx_orig <= mx_orig[i]+0.1))[0]
                         
    # check if duplicate
    if len(eqidx) > 1:
        delTrue = True
        
        # check if any authorities == EHB
        ehbidx = where(auth[eqidx] == 'EHB')[0]
        
        # just use first rec and delete rest
        if len(ehbidx) < 1:
            delidx = hstack((delidx, eqidx[1:]))
            	
        # else if ehb exists, use it
        else:
            notehb = where(auth[eqidx] != 'EHB')[0]
            delidx = hstack((delidx, eqidx[notehb]))
        
        # increment to next event
        i += len(eqidx) - 1
        
        print datestr[eqidx[0]], eqidx
    
    i += 1
    
delidx = array(delidx, dtype=int)

# make new lines
newlines = array(lines)

# delete duplicate lines - inc delidx because of header row
newlines = delete(newlines, delidx+1)

f = open(nsha2012csv.strip('.csv')+'.del.csv', 'wb')
f.writelines(newlines)
f.close()


