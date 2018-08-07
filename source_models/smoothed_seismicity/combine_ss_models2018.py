"""This code combines the smoothed seismicity models applying different b-values and mmax's in different tectonic regions
"""

import os, sys
import numpy as np
import ogr
import shapefile
from shapely.geometry import Point, Polygon

from source_models.logic_trees import logic_tree
from source_models.utils.pt2fault_distance import read_pt_source, combine_pt_sources
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.mfd.evenly_discretized import EvenlyDiscretizedMFD

from utilities import params_from_shp

def gr2inc_mmax(mfd, mmaxs, weights, model_weight=1.):
    """Function to convert a GR distribution to incremental MFD and 
    collapse Mmax logic tree branches
    :params model_weight:
        weight of total model, i.e. multiply all rates by this weight
    """
    mfd_type = type(mfd).__name__
  #  print mfd_type
    if mfd_type != 'TruncatedGRMFD':
        msg = 'Input MFD should be of type TruncatedGRMFD, found type %s' % mfd_type
        raise TypeError(msg)
    # Ensure we get rates for all mmax values
    mfd.max_mag = max(mmaxs)
    mag_bins, rates = zip(*mfd.get_annual_occurrence_rates())
    mag_bins = np.array(mag_bins)
    rates = np.array(rates)
   # print 'mag_bins'
   # print 'rates'
    new_rates = np.zeros(len(mag_bins))
    for mmax, weight in zip(mmaxs, weights):
        idx = np.where(np.isclose(mag_bins,(mmax-0.05), rtol=1e-2))[0][-1]+1
        new_rates[:idx] += rates[:idx]*weight
    new_rates = new_rates*model_weight
    new_mfd = EvenlyDiscretizedMFD(mfd.min_mag, mfd.bin_width, list(new_rates))
    return new_mfd
        

def combine_ss_models(filename_stem, domains_shp, params,lt, bval_key, output_dir='./',
                      nrml_version = '04', weight=1.):#, id_base = 'ASS'):
    """ Combine smoothed seismicity models based on tectonic region types
    :params filename_stem:
        String for the start of the xml filename for the source model,
        assuming generic components (non generic are inferred, 
        e.g. bvalue and completeness model)
    :params domains_shp:
        shapefile defining tectonic domain regions
    :params params:
        list of dicts containing parameters derivded from the shapefile
     :bval_key
         key for the dicts in params  as we are merging by 
         bvalues  (best, lower, upper)
    :params lt:
        LogicTree object containing relevant values and weights for Mmax
    :params outfile:
        output nrml formatted file
    """

    dsf = shapefile.Reader(domains_shp)
    dom_shapes = dsf.shapes()    
    # Get indicies of relevant fields
    for i, f in enumerate(dsf.fields):
        if f[0]=='CODE':
            code_index = i-1
        if f[0]=='TRT':
            trt_index = i-1

    hypo_depth_dist_nc = PMF([(0.5, 10.0),
                              (0.25, 5.0),
                              (0.25, 15.0)])
    hypo_depth_dist_c = PMF([(0.5, 5.0),
                             (0.25, 2.5),
                             (0.25, 10.0)])
    hypo_depth_dist_ex = hypo_depth_dist_c
    hypo_depth_dict = {'Cratonic': hypo_depth_dist_c,
                       'Non_cratonic': hypo_depth_dist_nc,
                       'Extended': hypo_depth_dist_ex}
    # FIXME! - Temporary solution until nodal plan logic tree
    # info can be read directly from shapefile attributes
    nodal_plane_dist = PMF([(0.3, NodalPlane(0, 30, 90)),
                            (0.2, NodalPlane(90, 30, 90)),
                            (0.3, NodalPlane(180, 30, 90)),
                            (0.2, NodalPlane(270, 30, 90))])

    merged_pts = []
 
    # Get mmax values and weights
    mmaxs = {}
    mmaxs_w = {}
    for dom in params:
        print 'Processing source %s' % dom['CODE']
        print dom['TRT']
        if dom['TRT'] == 'NCratonic':
            dom['TRT'] = 'Non_cratonic'
        # For the moment, only consider regions within AUstralia
        if dom['TRT'] == 'Active' or dom['TRT'] == 'Interface' or \
                dom['TRT'] == 'Oceanic' or \
                dom['TRT'] == 'Intraslab' or dom['CODE'] == 'NECS' or \
                dom['CODE'] == 'NWO': 
            print 'Source %s not on continental Australia, skipping' % dom['CODE']
            continue
        elif dom['TRT'] == 'Cratonic':
            if dom['DOMAIN'] == 1:
                mmax_values, mmax_weights = lt.get_weights('Mmax', 'Archean')
            else:
                mmax_values, mmax_weights = lt.get_weights('Mmax', 'Proterozoic')
#        elif dom['TRT'] == 'Active':
#            print 'MMax logic tree not yet defined for active crust, using extended crust'
#            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Extended')
        else:
            mmax_values, mmax_weights = lt.get_weights('Mmax', dom['TRT'])
        mmax_values = [float(i) for i in mmax_values]
        mmax_weights = [float(i) for i in mmax_weights]
        print mmax_values
        print mmax_weights
        mmaxs[dom['CODE']] = mmax_values
        mmaxs_w[dom['CODE']] = mmax_weights

        pt_ids = []
    #for trt, filename in filedict.iteritems():
    #    print trt
        completeness_string = 'comp'
        for ym in dom['COMPLETENESS']:
            completeness_string += '_%i_%.1f' % (ym[0], ym[1])
        mmin = dom['COMPLETENESS'][0][1]
        filename = "%s_b%.3f_mmin%.1f_%s.xml" % (
            filename_stem, dom[bval_key], mmin,
            completeness_string)
        print 'Parsing %s' % filename
        
        # TA kluge - hardwire jdg547 path
        jdgpath = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/'
        
        print filename
        
        # Only keep points within domain
        pts = read_pt_source(filename)
        
#        shapes = np.where(trt_types
        for shape in dsf.shapeRecords():
#            print code_index
            print shape.record[code_index]
            if shape.record[code_index] == dom['CODE']:
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
         #           if strikes[i] >=360:
         #               strikes[i]=strikes[i]-360
                nodal_plan_dist = PMF([(0.34, NodalPlane(strikes[0], 30, 90)),
                                       (0.34, NodalPlane(strikes[1], 30, 90)),
                                       (0.08, NodalPlane(strikes[2], 30, 90)),
                                       (0.08, NodalPlane(strikes[3], 30, 90)),
                                       (0.08, NodalPlane(strikes[4], 30, 90)),
                                       (0.08, NodalPlane(strikes[5], 30, 90))])
                if dom['CODE'] == 'WARM' or dom['CODE'] == 'WAPM':
                    print 'Define special case for WARM'
                    nodal_plan_dist = PMF([(0.75, NodalPlane(45, 90, 0)),
                                           (0.125, NodalPlane(strikes[0], 30, 90)),
                                           (0.125, NodalPlane(strikes[1], 30, 90))])
                if dom['CODE'] == 'FMLR':
                    print 'Define special case for FMLR, 0.5 thrust, 0.5 SS'
                    nodal_plan_dist = PMF([(0.17, NodalPlane(strikes[0], 30, 90)),
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
                dom_poly = Polygon(shape.shape.points)
                for pt in pts:
                    pt_loc = Point(pt.location.x, pt.location.y)
                    if pt_loc.within(dom_poly):
                        pt.tectonic_region_type = dom['TRT']
                        pt.nodal_plane_distribution = nodal_plane_dist # FIXME! update based on data extracted from shapefile
                        pt.hypocenter_distribution = hypo_depth_dist
                        pt.rupture_aspect_ratio=2
                        mfd = pt.mfd
                        new_mfd = gr2inc_mmax(mfd, mmaxs[dom['CODE']], mmaxs_w[dom['CODE']], weight)
                        pt.mfd = new_mfd
                        if pt.source_id in pt_ids:
                            print 'Point source %s already exists!' % pt.source_id
                            print 'Skipping this source for trt %s' % zone_trt
                        else:
                            merged_pts.append(pt)
                            pt_ids.append(pt.source_id)

    outfile = "%s_%s.xml" % (
            filename_stem, bval_key)
    outfile = os.path.join(output_dir, outfile)
    name = outfile.rstrip('.xml')
    if nrml_version == '04':
        nodes = list(map(obj_to_node, sorted(merged_pts)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        with open(outfile, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
    return outfile
    
            
if __name__ == "__main__":
#    filedict = {'Non_cratonic': 'source_model_adelaide_pts.xml'}
    output_dir = 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    domains_shp = '../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
    lt  = logic_tree.LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
    params = params_from_shp(domains_shp, trt_ignore=['Interface', 'Active', 'Oceanic', 'Intraslab'])
    filename_stem = 'Australia_Adaptive_K3'
    bestb_xml = combine_ss_models(filename_stem, domains_shp, params, lt, bval_key='BVAL_BEST',
                                  output_dir=output_dir, nrml_version = '04', 
                                  weight=0.5)
    upperb_xml = combine_ss_models(filename_stem, domains_shp, params, lt, bval_key='BVAL_UPPER',
                                   output_dir=output_dir, nrml_version = '04',
                                   weight=0.3)
    lowerb_xml = combine_ss_models(filename_stem, domains_shp, params, lt, bval_key='BVAL_LOWER',
                                   output_dir=output_dir, nrml_version = '04',
                                   weight=0.2)

    # combine all pt source models
    point_source_list = [bestb_xml, upperb_xml, lowerb_xml]
    filepath = os.path.join(output_dir, output_dir+'.xml')
    name = output_dir

    # read list of files
    pt_source_model_list =[]
    for point_source_model in point_source_list:
        print 'Reading %s' % point_source_model
        pt_model = read_pt_source(point_source_model)
        pt_source_model_list.append(pt_model)
    combine_pt_sources(pt_source_model_list, filepath, name , nrml_version='04',
                       id_location_flag = 'location')

