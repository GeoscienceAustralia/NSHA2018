from tools.oq_tools import return_annualised_haz_curves
from numpy import interp, exp, log, array, around
from os import getcwd

probs = array([0.02,0.01375,0.01,0.00445,0.002,0.0021,0.001,0.0005,0.000404,0.0002,0.0001])

if getcwd().startswith('/nas'):
    hazcurvefile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_curve-mean-SA(1.0)_1.csv'
    hazcurvefile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_curve-mean-SA(0.2)_1.csv'
    hazcurvefile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_curve-mean-PGA_1.csv'
    sitelistfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/nsha_cities.csv'
else:
    hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-SA(0.2)_1.csv'
    hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-PGA_1.csv'
    sitelistfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'

###############################################################################
# parse site file
###############################################################################

lines = open(sitelistfile).readlines() 
places = []
place_lat = []
place_lon = []

for line in lines:
    dat = line.strip().split(',')
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))
    places.append(dat[2])

# get hazard data
siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

###############################################################################
# match cities
###############################################################################
citylist = []
for sd in siteDict:
    for place, plon, plat in zip(places, place_lon, place_lat):
        if around(plon, decimals=2) == around(sd['lon'], decimals=2) \
           and around(plat, decimals=2) == around(sd['lat'], decimals=2):
               citylist.append(place)
           
           

interpTXT = 'LON,LAT,' + ','.join(('P'+str(x) for x in probs)) + ',LOCATION\n'

# loop through site dict
for site, city in zip(siteDict, citylist):
    interpHaz = exp(interp(log(probs[::-1]), log(site['poe_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    	
    # make text
    interpTXT += ','.join((str('%0.2f' % site['lon']), str('%0.2f' % site['lat']))) + ',' \
                 + ','.join((str(x) for x in interpHaz)) + ',' + city + '\n'
    
f = open('interp_haz_curves.csv', 'wb')
f.write(interpTXT)
f.close()