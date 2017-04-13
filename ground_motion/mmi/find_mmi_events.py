# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:57:07 2017

@author: u56903
"""
from os import path
from catalogue.parsers import parse_NSHA2012_catalogue
import datetime as dt

#########################################################################
# parse calalogue & convert to HMTK
#########################################################################

# Use 2012 NSHA catalogue
nsha2012csv = path.join('..', '..', 'catalogue', 'data', 'AUSTCAT.MW.V0.12.csv')
nsha_dict = parse_NSHA2012_catalogue(nsha2012csv)

# now parse mmi file
mmievfile = 'mmi_mags.csv'
lines = open(mmievfile).readlines()

eds = []
edt = []
emw = []
enr = []
epl = []

for line in lines:
    dat = line.strip().split(',')
    
    eds.append(dat[0])
    edt.append(dt.datetime.strptime(dat[0],'%Y%m%d%H%M'))
    emw.append(dat[1])
    enr.append(dat[2])
    epl.append(dat[3])

outtxt = 'DATETIME,MMI_MW,MX_ORIG,MX_REVML,PREF_MW,AUTH\n'    
for nd in nsha_dict:
    for ed, ds, em in zip(edt, eds, emw):
        if ed >= nd['datetime'] - dt.timedelta(minutes=1) \
            and ed <= nd['datetime'] + dt.timedelta(minutes=1):
                print ed
                outtxt += ','.join((ds, em, str(nd['mx_orig']), \
                                  str(nd['mx_rev_ml']), str(nd['prefmag']), \
                                  str(nd['auth']))) + '\n'
                                  
f = open('auscat_mmi_merge.csv', 'wb')
f.write(outtxt)
f.close()
                
    