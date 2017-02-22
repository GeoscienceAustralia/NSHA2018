import shapefile
from sys import argv
from shapely.geometry import Point, Polygon
try:
    from tools.nsha_tools import get_field_data
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

'''
Script to identify overlapping source zones in polygon shapefiles

Usage:
    
    python check_shape_overlap.py <path to shapefile> <ID field>
    
e.g:
    python check_shape_overlap.py ARUP/shapefiles/ARUP_source_model.shp sub_zone
'''    

###############################################################################
# parse shapefile
###############################################################################

shpfile = argv[1]
code_field = argv[2]

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
    
# get src code
codes = get_field_data(sf, code_field, 'str')

###############################################################################
# find poly points within other polygons
###############################################################################

# loop through zones 1st time
for code1, poly1 in zip(codes, shapes):
    
    # loop through zones 2nd time
    for code2, poly2 in zip(codes, shapes):
        
        # now loop through all points in poly2
        for vlon, vlat in poly2.points:
            point = Point(vlon, vlat)
            
            # chech if point in poly1
            if point.within(Polygon(poly1.points)) and code1 != code2:
                print code2, 'is overlapping', code1, ':', vlon, vlat
