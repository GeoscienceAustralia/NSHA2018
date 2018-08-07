"""Build NRML openquake input file for 2018 update to Hall et al 2007 model

Reference: Hall, L., Dimer, F., and Somerville, P. 2007. A Spatially Distributed Earthquake Source Model
for Australia.
Australian Earthquake Engineering Society 2007 Conference

Files received from Valentina Koschatzky, Risk Frontiers

Jonathan Griffin
Geoscience Australia
April 2018
"""

import os
import glob
import numpy as np 
# To build source model
import shapefile
from shapely.geometry import Polygon
import shapely.geometry
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import write, NAMESPACE
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
try:
    from tools.nsha_tools import get_field_data, get_shp_centroid
except:
    print 'Add PYTHONPATH to NSHA18 root directory'
from source_models.smoothed_seismicity.combine_ss_models import gr2inc_mmax
from source_models.smoothed_seismicity.utilities import params_from_shp
from source_models.logic_trees import logic_tree

source_data_filename = 'AB_Values.shp'
nrml_version='04'
msr = Leonard2014_SCR()
tom = PoissonTOM(50) 

source_data = shapefile.Reader(source_data_filename)
shapes = source_data.shapes()
a_values = get_field_data(source_data, 'aVal', 'float')
b_values = get_field_data(source_data, 'bVal', 'float')
ids = get_field_data(source_data, 'ID', 'float')

# Default values - real values should be based on neotectonic domains
min_mag = 4.5
max_mag = 7.2
depth = 10.0


###############################################################################
# get neotectonic domain number from centroid
###############################################################################
# load domains shp
domains_shapefile = '../../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
#os.path.join('..','..','zones','shapefiles','Domains','Domains_NSHA18_single_Mc.shp')
dsf = shapefile.Reader(domains_shapefile)
lt  = logic_tree.LogicTree('../../../shared/seismic_source_model_weights_rounded_p0.4.csv')
params = params_from_shp(domains_shapefile, trt_ignore=['Interface', 'Active', 'Oceanic', 'Intraslab'])
#print params
# get domains
#neo_doms  = get_field_data(dsf, 'DOMAIN', 'float')
#dom_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
#dom_trt  = get_field_data(dsf, 'GMM_TRT', 'str')
#dom_dep  = get_field_data(dsf, 'DEP_BEST', 'float')

# get domain polygons
#dom_shapes = dsf.shapes()

print a_values
# Loop over each cell (defined as a polygon/square, so we extract centroid)
# Get indicies of relevant fields
for i, f in enumerate(dsf.fields):
    if f[0]=='CODE':
        code_index = i-1
    if f[0]=='TRT':
        trt_index = i-1
merged_pts = []
pt_ids = []
# Get mmax values and weights
mmaxs = {}
mmaxs_w = {}
for dom in params:
    print 'Doing %s' % dom['CODE']
    if dom['TRT'] == 'NCratonic':
        dom['TRT'] = 'Non_cratonic'
    if dom['TRT'] == 'Cratonic':
        if dom['DOMAIN'] == 1:
            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Archean')
        else:
            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Proterozoic')
    else:
        mmax_values, mmax_weights = lt.get_weights('Mmax', dom['TRT'])
    mmax_values = [float(i) for i in mmax_values]
    mmax_weights = [float(i) for i in mmax_weights]    
    mmaxs[dom['CODE']] = mmax_values
    mmaxs_w[dom['CODE']] = mmax_weights
    for dom_shape in dsf.shapeRecords():
        if dom_shape.record[code_index] == dom['CODE']:
            # Check for undefined depths (-999 values)
            if dom['DEP_BEST'] < 0:
                print 'Setting best depth to 10 km'
                dom['DEP_BEST']=10
            if dom['DEP_UPPER'] < 0:
                print 'Setting upper depth to 5 km'
                dom['DEP_UPPER']=5
            if dom['DEP_LOWER'] < 0:
                print 'Setting lower depth to 15 km'
                dom['DEP_LOWER']=15
            hypo_depth_dist = PMF([(0.5, dom['DEP_BEST']),
                                   (0.25, dom['DEP_LOWER']),
                                   (0.25, dom['DEP_UPPER'])])
                # Define nodal planes as thrusts except for special cases
            str1 = dom['SHMAX'] + 90.
            str2 = dom['SHMAX'] + 270.
            str3 = dom['SHMAX'] + dom['SHMAX_SIG'] + 90.
            str4 = dom['SHMAX']+ dom['SHMAX_SIG'] + 270.
            str5 = dom['SHMAX'] - dom['SHMAX_SIG'] + 90.
            str6 = dom['SHMAX'] - dom['SHMAX_SIG'] + 270.
            strikes = [str1,str2,str3,str4,str5,str6]
            for i,strike in enumerate(strikes):
                if strike >=360:
                    strikes[i]=strike-360
            nodal_plane_dist = PMF([(0.34, NodalPlane(strikes[0], 30, 90)),
                                   (0.34, NodalPlane(strikes[1], 30, 90)),
                                   (0.08, NodalPlane(strikes[2], 30, 90)),
                                   (0.08, NodalPlane(strikes[3], 30, 90)),
                                   (0.08, NodalPlane(strikes[4], 30, 90)),
                                   (0.08, NodalPlane(strikes[5], 30, 90))])
            if dom['CODE'] == 'WARM' or dom['CODE'] == 'WAPM':
                print 'Define special case for WARM'
                nodal_plane_dist = PMF([(0.75, NodalPlane(45, 90, 0)),
                                       (0.125, NodalPlane(strikes[0], 30, 90)),
                                       (0.125, NodalPlane(strikes[1], 30, 90))])
            if dom['CODE'] == 'FMLR':
                print 'Define special case for FMLR, 0.5 thrust, 0.5 SS'
                nodal_plane_dist = PMF([(0.17, NodalPlane(strikes[0], 30, 90)),
                                       (0.17, NodalPlane(strikes[1], 30, 90)),
                                       (0.04, NodalPlane(strikes[2], 30, 90)),
                                       (0.04, NodalPlane(strikes[3], 30, 90)),
                                       (0.04, NodalPlane(strikes[4], 30, 90)),
                                       (0.04, NodalPlane(strikes[5], 30, 90)),
                                       (0.17, NodalPlane(strikes[0], 90, 0)),
                                       (0.17, NodalPlane(strikes[1], 90, 0)),
                                       (0.04, NodalPlane(strikes[2], 90, 0)),
                                       (0.04, NodalPlane(strikes[3], 90, 0)),
                                       (0.04, NodalPlane(strikes[4], 90, 0)),
                                       (0.04, NodalPlane(strikes[5], 90, 0))])        
   
            dom_poly = Polygon(dom_shape.shape.points)
            for i,shape in enumerate(shapes):
                centroid = get_shp_centroid(shape.points)
                shapely_pt = shapely.geometry.Point(centroid[0], centroid[1])
                #print centroid
                if shapely_pt.within(dom_poly):
                    pt = Point(centroid[0], centroid[1], depth) # Openquake geometry Point
                    tectonic_region_type = dom['GMM_TRT']
                    nodal_plane_distribution = nodal_plane_dist # FIXME! update based on data extracted from shapefile
                    hypocenter_distribution = hypo_depth_dist
                    rupture_aspect_ratio=2
                    mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, a_values[i], b_values[i])
                    new_mfd = gr2inc_mmax(mfd, mmaxs[dom['CODE']], mmaxs_w[dom['CODE']], model_weight=1.)
                    mfd = new_mfd
                    name = 'Hall_%i' % ids[i]
                    source_id = name
                    if source_id in pt_ids:
                        print 'Point source %s already exists!' % source_id
                        print 'Skipping this source for trt %s' % dom['TRT']
                    else:
                        pt_source = PointSource(source_id, name, dom['GMM_TRT'],
                                                mfd, 2, msr, 1.5,
                                                tom, 0.1, 20.0, pt,
                                                nodal_plane_dist, 
                                                hypo_depth_dist)
                        merged_pts.append(pt_source)
                        pt_ids.append(source_id)
outfile = 'Hall2007_2018.xml'
name = outfile.rstrip('.xml')
if nrml_version == '04':
    nodes = list(map(obj_to_node, sorted(merged_pts)))
    source_model = Node("sourceModel", {"name": name}, nodes=nodes)
    with open(outfile, 'wb') as f:
        nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
