"""Build NRML openquake input file for Hall et al 2007 model

Reference: Hall, L., Dimer, F., and Somerville, P. 2007. A Spatially Distributed Earthquake Source Model
for Australia.
Australian Earthquake Engineering Society 2007 Conference

Jonathan Griffin
Geoscience Australia
February 2016
"""

import os
import glob
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
    print 'Point PYTHONPATH to NSHA18 root directory'

# Default values - real values should be based on neotectonic domains
min_mag = 4.5
max_mag = 7.2
depth = 10.0
trt = 'Non_cratonic'

#############################################
# Parse original data (b_values from Hall et al 2007 paper
region_b_values = {'RF_NW-5.xyz' : 0.86,'RF_OZ-5.xyz' : 0.82,
                   'RF_S-5.xyz' : 0.84, 'RF_SE-5.xyz' : 0.82,
                   'RF_SW-5.xyz' : 0.70}
original_source_data_files = glob.glob(os.path.join('original_source_data', '*.xyz'))
for filename in original_source_data_files:
    print 'Reading data from %s' % filename
    filebase = filename.split('/')[-1]
    data = np.genfromtxt(filename, skip_header = 2)
    # Get dx and dy
    f_in = open(filename, 'r')
    header1 = f_in.readline()
    header2 = f_in.readline()
    header2 = header2.split()
    dx = float(header2[1])
    dy = float(header2[2])
    print dx, dy
    dxdy = dx*dy
    # a values in paper defined as annual rate of M > 5 per 100 km^2
    # i.e. they are NA5, not a0
    # Rates are defined per 100 km^2
    # convert to OpenQuake a value defintion
    # From Andreas: NA0 = NA5*10**(5*b)
    # Then we reduce to the grid cell size
    # N0 = NA0*dx*dy
    # a0 = log10(N0)
    b_val = region_b_values[filebase]
    try:
        lons = np.append(lons, data[:,0])
        lats = np.append(lats, data[:,1])
#        a_vals = np.append(a_vals, data[:,2]*np.power(10, b_val*5.0)/100.)
        a_vals = np.append(a_vals, np.log10(data[:,2]*np.power(10,5*b_val)*dxdy))
        b_vals = np.append(b_vals, np.ones(data[:,2].size)*b_val)
    except NameError:
        lons = data[:,0]
        lats = data[:,1]
#        a_vals = np.log10((data[:,2] + b_val*5.0)*dxdy)
#        a_vals = data[:,2]*np.power(10, b_val*5.0)/100.
#       a_vals = np.log10(data[:,2]*np.power(10,5*b_val)*dxdy)
        a_vals = np.log10(data[:,2]*np.power(10,5*b_val)*dxdy)
        b_vals = np.ones(data[:,2].size)*b_val
print a_vals

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

###############################################################################
# get TRT, depth, from Leonard08
###############################################################################
# load domains shp
lsf = shapefile.Reader(os.path.join('..','..','zones','Leonard2008','shapefiles','LEONARD08_NSHA18.shp'))

# get domains
ltrt  = get_field_data(lsf, 'TRT', 'str')
ldep  = get_field_data(lsf, 'DEP_BEST', 'float')
# get domain polygons
l08_shapes = lsf.shapes()

# Build sources
print 'Building point sources'
source_list = []
source_models = [] # Use later for logic tree
for j in range(len(lons)):
    identifier = 'H_' + str(j)
    name = 'Hall2007_' + str(j)
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
source_model = mtkSourceModel(identifier=0, name='Hall2007',
                              sources = source_list)
print 'Writing to NRML'
outbase = 'Hall2007'
source_model_filename = outbase + '_source_model.xml'
source_model.serialise_to_nrml(source_model_filename)
source_models.append(source_model)

######################################
# Now write the source model logic tree file
######################################
print 'Writing logic tree file'
newxml = '<?xml version="1.0" encoding="UTF-8"?>\n'
newxml += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
newxml += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
newxml += '    <logicTree logicTreeID="lt1">\n'
newxml += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n'
newxml += '            <logicTreeBranchSet uncertaintyType="sourceModel"\n' \
    '                                branchSetID="bs1">\n\n'

# make branches
for i, branch in enumerate(source_models):
    newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
    newxml += '                    <uncertaintyModel>'+source_model_filename+'</uncertaintyModel>\n'
    newxml += '                    <uncertaintyWeight>'+str(1)+'</uncertaintyWeight>\n'
    newxml += '                </logicTreeBranch>\n\n'
    
newxml += '            </logicTreeBranchSet>\n'
newxml += '        </logicTreeBranchingLevel>\n'
newxml += '    </logicTree>\n'
newxml += '</nrml>'
        
# write logic tree to file
outxml = outbase + '_source_model_logic_tree.xml'
f = open(outxml,'w')
f.write(newxml)
f.close()
