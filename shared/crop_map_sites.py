from shapely.geometry import Point, Polygon
import shapefile
from numpy import array

# parse site file



shpfile = '../source_models/zones/shapefiles/Other/hazard_crop.shp'
csvfile = '../source_models/zones/2018_mw/Domains_multi_mc/results_maps_PGA/hazard_map-mean_mapping.csv'

# parse file
print 'Reading csv file...'
lines = open(csvfile).readlines()[2:]

lon = []
lat = []

for line in lines:
    dat = line.split(',')
    lon.append(float(dat[0]))
    lat.append(float(dat[1]))

lon = array(lon)
lat = array(lat)

# down sample lo/la
dlon = lon[range(0, len(lon), 10)]
dlat = lat[range(0, len(lat), 10)]
    
#parese shapefile
print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
poly = Polygon(shapes[0].points)

print 'Looping thru sites...'
# loop through points
inlo = []
inla = []
outtxt = ''

for lo, la in zip(lon, lat):
    pt = Point(lo, la)
    
    # if in poly, or every tenth point
    if pt.within(poly):
       
       outtxt += ','.join((str(lo), str(la))) + '\n'
       
# now add downsampled outside
for lo, la in zip(dlon, dlat):
    pt = Point(lo, la)
    
    # if in poly, or every tenth point
    if pt.within(poly) == False:
       
       outtxt += ','.join((str(lo), str(la))) + '\n'


# write to file
f = open('nsha18_map_sites.csv', 'wb')
f.write(outtxt)
f.close()