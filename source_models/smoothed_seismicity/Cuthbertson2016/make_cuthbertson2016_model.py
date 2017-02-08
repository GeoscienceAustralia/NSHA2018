"""Build NRML openquake input file for Cuthbertson 2016 model

Reference: Cuthbertson, R. 2016. Automatic determination of seismicity
rates in Australia. Australian Earthquake Engineering 
Society 2016 Conference, Nov 25-27, Melbourne, Vic

Jonathan Griffin
Geoscience Australia
February 2016
"""

import os
import numpy as np 
# To build source model
import shapefile
from shapely.geometry import Polygon
import shapely.geometry
from hmtk.sources.source_model import mtkSourceModel
from hmtk.sources.point_source import mtkPointSource
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid
except:
    print 'Add PYTHONPATH to NSHA18 root directory'

# Default values - real values should be based on neotectonic domains
min_mag = 4.5
max_mag = 7.2
depth = 10.0
trt = 'Non_cratonic'

#############################################
# Parse original data
original_source_data = 'original_source_data/Aus_7_rates_and_area.txt'
data = np.genfromtxt(original_source_data, delimiter = ' ', 
                     skip_header = 1, dtype = ("|S24", float, float, int))
#print data
lons = []
lats = []
#print data['f0']
for location in data['f0']:
    location = location.split('_')
    lons.append(float(location[1]))
    lats.append(float(location[2]))
a_vals = data['f1']/4. # Divide by 4 as cells are overlapping
b_vals = data['f2']

###############################################################################
# get neotectonic domain number from centroid
###############################################################################
# load domains shp
dsf = shapefile.Reader(os.path.join('..','..','zones','Domains','shapefiles','DOMAINS_NSHA18.shp'))
# get domains
neo_doms  = get_field_data(dsf, 'DOMAIN', 'float')
dom_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
# get domain polygons
dom_shapes = dsf.shapes()
#n_dom = []
#n_mmax = []

###############################################################################
# get TRT, depth from Leonard08
###############################################################################
# load domains shp
lsf = shapefile.Reader(os.path.join('..','..','zones','Leonard2008','shapefiles','LEONARD08_NSHA18.shp'))

# get domains
ltrt  = get_field_data(lsf, 'TRT', 'str')
ldep  = get_field_data(lsf, 'DEP_BEST', 'float')
# get domain polygons
l08_shapes = lsf.shapes()

# Build sources
source_list = []
for j in range(len(lons)):
    identifier = 'RC_' + str(j)
    name = 'Cuthbertson_' + str(j)
    # Need to use shapely functions intially
    shapely_pt = shapely.geometry.Point(lons[j], lats[j])
    # Get parameters based on domain
    # loop through domains and find point in poly
    for neo_dom, mmax, dom_shape in zip(neo_doms, dom_mmax, dom_shapes):
        dom_poly = Polygon(dom_shape.points)       
        # check if leonard centroid in domains poly
        if shapely_pt.within(dom_poly):
            tmp_dom = neo_dom
            max_mag = mmax
    for zone_trt, zone_dep, l_shape in zip(ltrt, ldep, l08_shapes):
        l_poly = Polygon(l_shape.points)
        if shapely_pt.within(l_poly):
            trt = zone_trt
            depth = zone_dep
#    print max_mag, trt, depth                       

    point = Point(lons[j], lats[j], depth) # Openquake geometry Point
    mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, a_vals[j], b_vals[j])
    hypo_depth_dist = PMF([(1.0, depth)])
    nodal_plane_dist = PMF([(0.3, NodalPlane(0, 30, 90)),
                            (0.2, NodalPlane(90, 30, 90)),
                            (0.3, NodalPlane(180, 30, 90)),
                            (0.2, NodalPlane(270, 30, 90))])
    point_source = mtkPointSource(identifier, name, geometry=point, mfd=mfd,
                           mag_scale_rel = 'Leonard2014_SCR', rupt_aspect_ratio=1.0,
                           upper_depth = 0.1, lower_depth = 20.0,
                           trt = trt, nodal_plane_dist = nodal_plane_dist,
                           hypo_depth_dist = hypo_depth_dist)
    source_list.append(point_source)
source_model = mtkSourceModel(identifier=0, name='Cuthbertson2016',
                              sources = source_list)
source_model.serialise_to_nrml('source_model_cuthbertson2016.xml')
