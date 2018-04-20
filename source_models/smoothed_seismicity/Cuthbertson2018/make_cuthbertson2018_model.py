"""Build NRML openquake input file for Cuthbertson 2016 model

Reference: Cuthbertson, R. 2016. Automatic determination of seismicity
rates in Australia. Australian Earthquake Engineering 
Society 2016 Conference, Nov 25-27, Melbourne, Vic

Jonathan Griffin
Geoscience Australia
February 2016
"""

import os, sys
import glob
import numpy as np 
from scipy import interpolate
# To build source model
import shapefile
from shapely.geometry import Polygon
import shapely.geometry
from hmtk.sources.source_model import mtkSourceModel
from hmtk.sources.point_source import mtkPointSource
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

# Default values - real values should be based on neotectonic domains
min_mag = 4.5
max_mag = 7.2
depth = 10.0
trt = 'Non_cratonic'

#############################################
# Parse original data
original_source_data = 'Aus_9_rates_and_area.txt'
print 'Reading data from %s' % original_source_data
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
rates = data['f1']/4.
a_vals = np.log10(data['f1']/4.) # Divide by 4 as cells are overlapping
b_vals = data['f2']

###############################################################################
# get neotectonic domain number from centroid
###############################################################################
# load domains shp
domains_shp =  '../../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
#os.path.join('..','..','zones','shapefiles','Domains','Domains_NSHA18_single_Mc.shp')
dsf = shapefile.Reader(domains_shp)
lt  = logic_tree.LogicTree('../../../shared/seismic_source_model_weights_rounded_p0.4.csv')
params = params_from_shp(domains_shp, trt_ignore=['Interface', 'Active', 'Oceanic', 'Intraslab'])
# get domains
neo_doms  = get_field_data(dsf, 'DOMAIN', 'float')
dom_codes = get_field_data(dsf, 'CODE', 'str')
dom_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
dom_trt  = get_field_data(dsf, 'GMM_TRT', 'str')
dom_dep  = get_field_data(dsf, 'DEP_BEST', 'float')

# Build some dictionaries of parameters for each domain
param_index = {}
k = 0
print params
for dom in params:
    print 'Processing source %s' % dom['CODE']
    if dom['TRT'] == 'NCratonic' or dom['TRT'] == 'Extended':
        dom['TRT'] = 'Non_cratonic'
        # For the moment, only consider regions within AUstralia
    if dom['TRT'] == 'Active' or dom['TRT'] == 'Interface' or \
            dom['TRT'] == 'Oceanic' or \
            dom['TRT'] == 'Intraslab' or dom['CODE'] == 'NECS' or \
            dom['CODE'] == 'NWO': 
        print 'Source %s not on continental Australia, skipping' % dom['CODE']
        k+=1
        continue
    elif dom['TRT'] == 'Cratonic':
        if dom['DOMAIN'] == 1:
            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Archean')
        else:
            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Proterozoic')
    else:
        mmax_values, mmax_weights = lt.get_weights('Mmax', dom['TRT'])
    mmax_values = [float(i) for i in mmax_values]
    mmax_weights = [float(i) for i in mmax_weights]
    dom['MMAXS'] = mmax_values
    dom['MMAX_WEIGHTS'] = mmax_weights

    param_index[dom['CODE']]=k
    k+=1  

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
            #           if strikes[i] >=360:
            #               strikes[i]=strikes[i]-360
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
    dom['DEPTH_DIST'] = hypo_depth_dist
    dom['NP_DIST'] = nodal_plane_dist
# get domain polygons
dom_shapes = dsf.shapes()

###############################################################################
# get TRT, depth from Leonard08
###############################################################################

# Build sources
print 'Building point sources'
source_list = []
source_models = [] # Use later for logic tree
degoff = np.arange(-0.2,0.3,0.1) # loc offset for distributing rates

# Interpolation of a values
#rates_split = a_vals /  25. # Normalisation
#rates_split = rates/25.
#new_lons = []
#new_lats = []
#a_split=[]
#for lon,lat,rate in zip(lons,lats,rates):
#    for dlo in degoff:
#        for dla in degoff:
#            new_lons.append(lon+dlo)
#            new_lats.append(lat+dla)
           # a_split.append(np.log10(rate/25.))
#print len(lons)
#print len(lats)
#print len(rates_split)
#a_split = np.log10((10**a_vals) /  25.) # Normalisation
##print a_split
#interp_fn = interpolate.interp2d(lons, lats, a_split, kind='linear', fill_value=1.)
#interp_rates = interp_fn(new_lons,new_lats)
#print interp_rates
#sys.exit
#a_split = np.log10((10**a_vals) /  25.) # Normalisation


print param_index
for j in range(len(lons)):
    identifier = 'RC_' + str(j)
    name = 'Cuthbertson_' + str(j)
    # Need to use shapely functions intially
    shapely_pt = shapely.geometry.Point(lons[j], lats[j])
    # Get parameters based on domain
    # loop through domains and find point in poly
    for neo_dom, dom_code, zone_trt, zone_dep, dom_shape in zip(neo_doms, dom_codes, dom_trt, dom_dep, dom_shapes):
        dom_poly = Polygon(dom_shape.points)       
        # check if leonard centroid in domains poly
        if shapely_pt.within(dom_poly):
            tmp_dom = dom_code
#            print tmp_dom
#            print zone_trt
            try:
                k = param_index[tmp_dom]
            except KeyError:
                print 'Skipping points in offshore zone %s' % tmp_dom
                cont = True
                continue
#            max_mag = mmax
            trt = zone_trt
#            depth = zone_dep
    if cont:
        cont=False
        continue
    # split rates 25 ways
    a_split = np.log10((10**a_vals[j]) /  25.)
    
    # now loop through subdivided cells
    inc = 0
    for dlo in degoff:
        for dla in degoff:
            inc += 1
            splitIdentifier = identifier+'-'+str(inc)
            
            # set lola offset
            lo = lons[j]+dlo
            la = lats[j]+dla
            
            point = Point(lo, la, depth) # Openquake geometry Point
            
            mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, a_split, b_vals[j])
            weight = 1.0 # We aren't combining models here 
#            print params[k]['MMAXS']
#            print params[k]['CODE']
            new_mfd = gr2inc_mmax(mfd, params[k]['MMAXS'], params[k]['MMAX_WEIGHTS'], weight)
            
            hypo_depth_dist = params[k]['DEPTH_DIST'] #PMF([(1.0, depth)])
#            print hypo_depth_dist
#            nodal_plane_dist = PMF([(0.3, NodalPlane(0, 30, 90)),
#                                    (0.2, NodalPlane(90, 30, 90)),
#                                    (0.3, NodalPlane(180, 30, 90)),
#                                    (0.2, NodalPlane(270, 30, 90))])
            nodal_plane_dist =  params[k]['NP_DIST']
#            print nodal_plane_dist
            point_source = mtkPointSource(splitIdentifier, name, geometry=point, mfd=new_mfd,
                                   mag_scale_rel = 'Leonard2014_SCR', rupt_aspect_ratio=1.5,
                                   upper_depth = 0.1, lower_depth = 20.0,
                                   trt = trt, nodal_plane_dist = nodal_plane_dist,
                                   hypo_depth_dist = hypo_depth_dist)
            
            source_list.append(point_source)

source_model = mtkSourceModel(identifier=0, name='Cuthbertson2018',
                              sources = source_list)

print 'Writing to NRML'
outbase = 'cuthbertson2018'
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
