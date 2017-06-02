# -*- coding: utf-8 -*-
"""
Created on Thu Apr 06 09:26:34 2017

@author: u56903
"""
import datetime as dt
from numpy import nan

eafile = 'EA_earthquakes_with_MW.txt'
wafile = 'F:\Catalogues\ML2MW\WA_earthquakes_with_MW.txt'
hgfile = 'M:\NSHA2018\catalogue\data\kiwi_results.dat'

evdt = []
evlo = []
evla = []
evdp = []
evmw = []
msrc = []

# read ea file
lines = open(eafile).readlines()[1:]

for line in lines:
    dat = line.strip().split('\t')
    sdt = dat[0].split()
    if len(sdt[1]) == 7:
        sdt[1] = '0'+sdt[1]

    # now get dt
    tdt = dt.datetime.strptime(' '.join(sdt), '%Y-%m-%d %H:%M:%S') 
    
    # append data
    evdt.append(tdt)
    evlo.append(float(dat[1]))
    evla.append(float(dat[2]))
    evdp.append(float(dat[3]))
    evmw.append(float(dat[4]))
    msrc.append(dat[5])
              
# read wa file
lines = open(wafile).readlines()[1:]

for line in lines:
    dat = line.strip().split('\t')

    # now get dt
    tdt = dt.datetime.strptime(dat[0], '%Y-%m-%d %H:%M:%S') 
    
    # append data
    evdt.append(tdt)
    evlo.append(float(dat[1]))
    evla.append(float(dat[2]))
    evdp.append(float(dat[3]))
    evmw.append(float(dat[4]))
    msrc.append(dat[6])

# read hg file
lines = open(hgfile).readlines()[1:]

for line in lines:
    dat = line.strip().split(' ')
    
    dta = [int(round(float(x))) for x in dat[3:9]]
    # now get dt
    tdt = dt.datetime(dta[0], dta[1], dta[2], dta[3], dta[4], dta[5]) 
    
    # append data
    evdt.append(tdt)
    evlo.append(float(dat[2]))
    evla.append(float(dat[1]))
    evdp.append(nan)
    evmw.append(float(dat[9]))
    msrc.append('Ghasemi et al (2016)')

# Add 1941 Meeberrie

# get dt
tdt = dt.datetime(1941, 4, 29, 1, 35, 38.57) 

# append data
evdt.append(tdt)
evlo.append(116.248)
evla.append(-26.826)
evdp.append(15.4)
evmw.append(6.18)
msrc.append('GEM-ISC')

# write to file
outtxt = 'DATETIME,LON,LAT,DEP,MW,SRC\n'

for i in range(0,len(evdt)):
    newline = ','.join((evdt[i].strftime('%Y-%m-%d %H:%M:%S'), str('%0.3f' % evlo[i]), \
                        str('%0.3f' % evla[i]), str('%0.1f' % evdp[i]), \
                        str('%0.2f' % evmw[i]), msrc[i])) + '\n'
    
    outtxt += newline
                        
f = open('combined_au_mw.dat', 'wb')
f.write(outtxt)
f.close()
