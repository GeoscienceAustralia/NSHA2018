from tools.oq_tools import return_annualised_haz_curves
from numpy import interp, exp, log, array, around
from os import getcwd, path, mkdir
from sys import argv

probs = array([0.02,0.01375,0.01,0.00445,0.002,0.0021,0.001,0.0005,0.000404,0.0002,0.0001])

periods = ['PGA', 'SA(0.05)', 'SA(0.07)', 'SA(0.1)', 'SA(0.2)', 'SA(0.3)', 'SA(0.4)', \
           'SA(0.1)']

hazcurvefile = argv[1]

# check to see if maps exists
if path.isdir('hazard_grids') == False:
    mkdir('hazard_grids')

###############################################################################
# parse hazard data
###############################################################################

# get hazard data
siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)           

interpTXT = 'LON,LAT,' + ','.join(('P'+str(x) for x in probs)) + '\n'

# loop through site dict
for site in siteDict:
    interpHaz = exp(interp(log(probs[::-1]), log(site['poe_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    	
    # make text
    interpTXT += ','.join((str('%0.2f' % site['lon']), str('%0.2f' % site['lat']))) + ',' \
                 + ','.join((str(x) for x in interpHaz)) + '\n'

###############################################################################
# make output
###############################################################################

outfile = 

outpath = path.join('hazard_grids', 'interp_haz_curves.csv')
f = open(outcsv, 'wb')
f.write(interpTXT)
f.close()