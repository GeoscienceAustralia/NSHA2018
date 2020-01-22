# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 12:25:12 2018

@author: u56903
"""
from tools.oq_tools import return_annualised_haz_curves
from numpy import array, mgrid, linspace, isinf, hstack, interp, log, exp
from misc_tools import dictlist2array
from scipy.interpolate import griddata
from os import path, mkdir
from hazard_tools import get_percent_chance_from_return_period, \
                         get_probability_from_percent_chance

##############################################################################
# set some default values here
##############################################################################

interp_lon = array([116.17])
interp_lat = array([-31.90])

hazcurvefile = 'results_maps_PGA/hazard_curve-mean-PGA_1.csv'

##############################################################################
# parse hazard grid
##############################################################################

# parse grid file
gridDict, poe_imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

lons = dictlist2array(gridDict, 'lon')
lons = lons.reshape(len(lons), 1)

lats = dictlist2array(gridDict, 'lat')
lats = lats.reshape(len(lats), 1)

localities = hstack((lons, lats))

##############################################################################
# interpolate hazard grid for each iml
##############################################################################

# for each IML, grid and interpolate
interp_poe = []

hazcurvetxt = 'IML, Annual PoE\n'
for i, iml in enumerate(poe_imls):
    # for each iml, get PoE
    poe = []
    for gd in gridDict:
        if isinf(gd['poe_probs_annual'][i]):
            # set to very small number
            poe.append(1E-20)
        else:
            poe.append(gd['poe_probs_annual'][i])
            #poe.append(gd['poe_probs_invtime'][i])
        
    poe = array(poe)
    
    # grid the data.
    interp_poe.append(griddata(localities, poe, (interp_lon, interp_lat), method='cubic')[0])
    hazcurvetxt += ','.join((str(iml), str('%0.5e' % interp_poe[-1]))) + '\n'

interp_poe = array(interp_poe)
##############################################################################
# export hazard curve
##############################################################################
    
# check to see if maps exists
if path.isdir('interp_hazard') == False:
    mkdir('interp_hazard')

# set filename
outhazcurve = '_'.join(('hazard_curve-mean-' + hazcurvefile.split('-')[-1].split('_')[0], \
                        str(interp_lon[0]), str(interp_lat[0])+'.csv'))

# write to file
f = open(path.join('interp_hazard', outhazcurve), 'wb')
f.write(hazcurvetxt)
f.close()

##############################################################################
# interpolate to standard return periods
##############################################################################

# set grid return periods for hazard curve
return_periods = ['100', '250', '475', '500', '800', '1000', '1500', '2000', '2475', \
                  '2500', '3000', '5000']


haztxt = 'RETURN_PERIOD,ANNUAL_PROBABILITY,HAZARD_LEVEL(g)\n'
for return_period in return_periods:
    percent_chance = get_percent_chance_from_return_period(float(return_period), investigation_time)
    return_period_num, probability = get_probability_from_percent_chance(percent_chance, investigation_time)
    
    interphaz = exp(interp(log(probability), log(interp_poe[::-1]), log(poe_imls[::-1])))
    haztxt += ','.join((return_period, str('%0.4e' % probability), str('%0.4e' % interphaz))) + '\n'


# set filename
outhazcurve = '_'.join(('return-period-mean-' + hazcurvefile.split('-')[-1].split('_')[0], \
                        str(interp_lon[0]), str(interp_lat[0])+'.csv'))

# write to file
f = open(path.join('interp_hazard', outhazcurve), 'wb')
f.write(haztxt)
f.close()




















    
