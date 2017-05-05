from __future__ import unicode_literals
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
from numpy import array, exp, log, interp, vstack, mean, around
from tools.oq_tools import return_annualised_haz_curves
from os import path, mkdir
import warnings, sys
reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

###############################################################################
# read config file
###############################################################################

conf_file = sys.argv[1]

# get paths for input files
lines = open(conf_file).readlines()
hazcurvefile1 = lines[0].split('=')[-1].strip()
job1 = hazcurvefile1.split('/')[0]
hazcurvefile2 = lines[1].split('=')[-1].strip()
job2 = hazcurvefile2.split('/')[0]
outputdir     = lines[2].split('=')[-1].strip()
sitelistfile  = lines[3].split('=')[-1].strip()
period        = lines[4].split('=')[-1].strip()

# check to see if exists
if path.isdir(outputdir) == False:
    mkdir(outputdir)

###############################################################################
# read OQ data
###############################################################################

def get_oq_haz_curves(hazcurvefile):
    if hazcurvefile.endswith('xml'):
        # Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
        lines = open(hazcurvefile, 'r').readlines()
        lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
        out = open(hazcurvefile, 'w')
        out.writelines(lines)
        out.close()
    
    # get annualize the curves.
    curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
    imls = array(metadata['imls'])
    
    return curves, curlon, curlat, metadata, imls

# get curves
curves1, curlon1, curlat1, metadata1, imls1 = get_oq_haz_curves(hazcurvefile1)
curves2, curlon2, curlat2, metadata2, imls2 = get_oq_haz_curves(hazcurvefile2)

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

###############################################################################
# write OQ & Frisk data to file
###############################################################################
i = 1
ii = 0
yhaz = 1./2475.

# set NBCC probabilities
nbccprobs = array([0.02,  0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001])

# set headers
if period == 'PGA' or period == 'PGV':
    jobhead = ' '.join(('Mean', period, 'Hazard (g)'))
else:
    jobhead = ' '.join(('Mean','SA['+period+']', 'Hazard (g)'))
    
oqhead1 = '\n \n--- '+job1+' ---\n \nLocation,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
oqhead2 = '\n \n--- '+job2+' ---\n \nLocation,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
rathead = '\n \n--- Ratios --- \t'+job1+' / '+job2+'\n \nLocation,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
pcdhead = '\n \n--- % Difference --- \t'+job1+' / '+job2+'\n \nLocation,P 0.02,P 0.01375,P 0.01,P 0.00445,P 0.0021,P 0.001,P 0.0005,P 0.000404,P 0.0002,P 0.0001\n'
oqt1 = ''
oqt2  = ''
rat = ''
pcd = ''

# loop thru OQ curves
# loop thru 1st OQ curves
for lon1, lat1, curve1 in zip(curlon1, curlat1, curves1):
    
    # loop thru 2nd OQ curves
    for lon2, lat2, curve2 in zip(curlon2, curlat2, curves2):
        if around(lon2, decimals=2) == around(lon1, decimals=2) \
           and around(lat2, decimals=2) == around(lat1, decimals=2):
            
            # interp to nbcc probs
            oqinterp1 = exp(interp(log(nbccprobs[::-1]), log(curve1[::-1]), log(imls1[::-1])))[::-1]
            oqinterp2 = exp(interp(log(nbccprobs[::-1]), log(curve2[::-1]), log(imls2[::-1])))[::-1]
            
            # loop thru places to get title
            for place, plon, plat in zip(places, place_lon, place_lat):
                if around(plon, decimals=2) == around(lon1, decimals=2) \
                   and around(plat, decimals=2) == around(lat1, decimals=2):
                    wrtplace = place
                    
            # set OQ1 text
            oqt1 += wrtplace + ',' + ','.join((str('%0.4f' % x) for x in oqinterp1)) + '\n'
                    
            # set frisk test
            oqt2 += wrtplace + ',' + ','.join((str('%0.4f' % x) for x in oqinterp2)) + '\n'
            
            # get ratios
            hazrat = oqinterp1 / oqinterp2
            
            # set ratio text
            rat += wrtplace + ',' + ','.join((str('%0.4f' % x) for x in hazrat)) + '\n'
            
            # calc % difference
            numer = oqinterp2 - oqinterp1
            denom = mean(vstack((oqinterp1, oqinterp2)), axis=0)
            pcdiff = 100. * (numer / denom)
            
            # set % diff text
            pcd += wrtplace + ',' + ','.join((str('%0.2f' % x) for x in pcdiff)) + '\n'
            

# combine text
outtxt = jobhead + oqhead1 + oqt1 + oqhead2 + oqt2 + rathead + rat + pcdhead + pcd

# write to file
csvfile =  path.join(outputdir, job1 +'_'+ job2 + '_hazard_ratio_' + str(period) +'.csv')
f = open(csvfile, 'wb')
f.write(outtxt)
f.close()
            
            
            
            
            
            
            
            
            
            
            
            
