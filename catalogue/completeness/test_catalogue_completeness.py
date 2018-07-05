from tools.mfd_tools import get_events_in_poly, parse_hmtk_cat
from tools.mapping_tools import get_field_data
from os import path
import shapefile
from shapely.geometry import Point, Polygon
from numpy import arange, logspace, nan, vstack, where
import matplotlib.pyplot as plt

################################################################################
# set params
depmin = -999
depmax = 999
mrng = arange(2.5, 6.6, 0.5)
trng = logspace(-1,3,119)[0:-20] # remove non-germane ranges
	
# set figures
ncolours = 10
cmap = plt.get_cmap('gist_ncar_r', ncolours)

# map colours
cs = (cmap(arange(ncolours)))
cs = cs[1:]

################################################################################
# first parse catalogue
hmtk_csv = path.join('..','data','NSHA18CAT_V0.1_hmtk_declustered.csv')
ggcat, neq = parse_hmtk_cat(hmtk_csv)

################################################################################
# now parse completeness zone shapefile
print 'Reading source shapefile...'
shpfile = path.join('completeness_zones.shp')
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
# get input arrays from shapefile
cmp_reg = get_field_data(sf, 'NAME', 'str')

################################################################################
# loop thru zones and extract earthquakes
n = 0
for poly, reg in zip(polygons, cmp_reg):
    n += 1
    plt.figure(n, figsize=(8,8))
    
    # get earthquakes in poly
    mvect, mxvect, tvect, dec_tvect, ev_dict \
                   = get_events_in_poly(ggcat, poly, depmin, depmax)
                   
    time2present = dec_tvect[-1] - dec_tvect
    
    # loop through mag ranges
    for j, m in enumerate(mrng):
        # loop through time ranges
        tmpRates = []
        for t in trng:
            # find events withing m & t ranges
            idx = where((mvect >= m) & (time2present <= t))[0] #& (mvect < m+0.5)
            
            if len(idx) == 0:
                tmpRates.append(nan)
            else:
                tmpRates.append(len(idx) / t)
                
        if m == mrng[0]:
            rates = tmpRates
        else:
            rates = vstack((rates, tmpRates))
            
        # plot completeness
        #label = ' '.join(('MW',str(m)+'-'+str(m+0.5)))
        label = ' '.join(('MW >=',str(m)))
        plt.loglog(trng, tmpRates, '-', c=cs[j], lw=2.0, label=label)
    
    # make pretty
    plt.grid(which='both')
    plt.xlabel('Time to Present (years)')
    plt.ylabel('Annual Number of Earthquakes (Mw >= Mx)')
    plt.title('Zone: '+reg)
    pltName = path.join('plots', reg+'.comp.png')
    plt.legend(loc=3, fontsize=13)
    plt.xlim([0.1, 300])
    plt.savefig(pltName, fmt='png', bbox_inches='tight')
    plt.close()
    