"""This code combines the smoothed seismicity models applying different b-values and mmax's in different tectonic regions
"""

import os, sys
import numpy as np
import ogr
import shapefile
from shapely.geometry import Point, Polygon

from NSHA2018.source_models.logic_trees import logic_tree
from NSHA2018.source_models.utils.pt2fault_distance import read_pt_source
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.mfd.evenly_discretized import EvenlyDiscretizedMFD

def gr2inc_mmax(mfd, mmaxs, weights):
    """Function to convert a GR distribution to incremental MFD and 
    collapse Mmax logic tree branches
    """
    mfd_type = type(mfd).__name__
    if mfd_type != 'TruncatedGRMFD':
        msg = 'Input MFD should be of type TruncatedGRMFD'
        raise(msg)
    # Ensure we get rates for all mmax values
    mfd.max_mag = max(mmaxs)
    mag_bins, rates = zip(*mfd.get_annual_occurrence_rates())
    mag_bins = np.array(mag_bins)
    rates = np.array(rates)
    new_rates = np.zeros(len(mag_bins))
    for mmax, weight in zip(mmaxs, weights):
        idx = np.where(np.isclose(mag_bins,(mmax-0.05), rtol=1e-2))[0][-1]+1
        new_rates[:idx] += rates[:idx]*weight
    new_mfd = EvenlyDiscretizedMFD(mfd.min_mag, mfd.bin_width, list(new_rates))
    return new_mfd
        

def combine_ss_models(filedict, domains_shp, lt, outfile,
                      nrml_version = '04'):
    """ Combine smoothed seismicity models based on tectonic region types
    :params filedict:
        dict of form filedict[trt] = filename specifying input file for that region
    :params domains_shp:
        shapefile defining tectonic domain regions
    :params lt:
        LogicTree object containing relevant values and weights for Mmax
    :params outfile:
        output nrml formatted file
    """
    print 'Getting tectonic region type from %s' % domains_shp
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(domains_shp, 0)
    dsf = data_source.GetLayer()
    trt_types = []
    for feature in dsf:
        trt_types.append(feature.GetField('TRT'))
    dsf = shapefile.Reader(domains_shp)
    dom_shapes = dsf.shapes()

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
    nodal_plane_dist = PMF([(0.3, NodalPlane(0, 30, 90)),
                            (0.2, NodalPlane(90, 30, 90)),
                            (0.3, NodalPlane(180, 30, 90)),
                            (0.2, NodalPlane(270, 30, 90))])

    merged_pts = []
    
    # Get mmax values and weights
    mmaxs = {}
    mmaxs_w = {}
    for trt, filename in filedict.iteritems():
        if trt == 'Cratonic':
            mmax_values, mmax_weights = lt.get_weights('Mmax', 'Proterozoic')
        else:
            mmax_values, mmax_weights = lt.get_weights('Mmax', trt)
        mmax_values = [float(i) for i in mmax_values]
        mmax_weights = [float(i) for i in mmax_weights]
        print mmax_values
        print mmax_weights
        mmaxs[trt] = mmax_values
        mmaxs_w[trt] = mmax_weights

    for trt, filename in filedict.iteritems():
        print trt
        print 'Parsing %s' % filename
        # Only keep points within domain
        pts = read_pt_source(filename)
#        shapes = np.where(trt_types
        for zone_trt, dom_shape in zip(trt_types, dom_shapes):
            dom_poly = Polygon(dom_shape.points)
            for pt in pts:
                pt_loc = Point(pt.location.x, pt.location.y)
                if pt_loc.within(dom_poly):
                    pt.tectonic_region_type = zone_trt
                    pt.nodal_plane_distribution = nodal_plane_dist
                    pt.hypocenter_distribution = hypo_depth_dict[zone_trt]
                    pt.rupture_aspect_ratio=2
                    mfd = pt.mfd
                    new_mfd = gr2inc_mmax(mfd, mmaxs[trt], mmaxs_w[trt])
                    pt.mfd = new_mfd
                    merged_pts.append(pt)
    
    name = outfile.rstrip('.xml')
    if nrml_version == '04':
        nodes = list(map(obj_to_node, sorted(merged_pts)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        with open(outfile, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)


if __name__ == "__main__":
#    filedict = {'Non_cratonic': 'source_model_adelaide_pts.xml'}
    filedict = {'Cratonic': 'source_model_Australia_Adaptive_K3_b0.819.csv.xml',
                'Non_cratonic': 'source_model_Australia_Adaptive_K3_b1.208.csv.xml',
                'Extended': 'source_model_Australia_Adaptive_K3_b0.835.csv.xml'}
    domains_shp = '../zones/2012_mw_ge_4.0/NSHA13_Background/shapefiles/NSHA13_BACKGROUND_NSHA18_MFD.shp'
    outfile = 'source_model_Australia_Adaptive_K3_merged.xml'
    lt  = logic_tree.LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
    combine_ss_models(filedict, domains_shp, lt, outfile, nrml_version = '04')
