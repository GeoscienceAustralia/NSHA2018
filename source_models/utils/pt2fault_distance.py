"""Reads an input point source model and generates ruptures according 
to the MFD and nodal plane distirbution. Nearest distance to the fault is calculated and used to reduce the Mmax of the point source model. 
"""

import os, sys
import numpy as np
import copy
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.sourceconverter import SourceConverter, SourceGroup
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
from openquake.hazardlib.geo.geodetic import distance, geodetic_distance
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
#plotting
from matplotlib import pyplot
#from hmtk.plotting.mapping import HMTKBaseMap


def read_pt_source(pt_source_file):
    """Read nrml source model into pt source objects
    """
    converter = SourceConverter(50, 2, width_of_mfd_bin=0.1,
                                area_source_discretization=10.)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(pt_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(pt_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
#    for pt in sources:
#        print pt.mfd.max_mag
    return sources

#def merge_pt_sources(point_source_list, filename, name, nrml_version='04'):
#    """Method for merging multiple point source file into one
#    """

def merge_rates(pt, added_pt, method='Add'):
    """Merge the rates of multiple mfds
    """
    if method == 'Add':
        mfd_type = type(pt.mfd).__name__
        if mfd_type == 'EvenlyDiscretizedMFD':
            mag_bins, rates = zip(*pt.mfd.get_annual_occurrence_rates())
            mag_bins = np.array(mag_bins)
            rates = np.array(rates)
            added_pt_mag_bins, added_pt_rates = zip(*added_pt.mfd.get_annual_occurrence_rates()) 
            added_pt_mag_bins = np.array(added_pt_mag_bins)
            added_pt_rates = np.array(added_pt_rates)
            new_rates = []
            for mag_bin in mag_bins:
                total_rate = np.sum(rates[np.where(mag_bins == mag_bin)]) + \
                    np.sum(added_pt_rates[np.where(added_pt_mag_bins == mag_bin)])
                new_rates.append(total_rate)
            return new_rates
    else:
        msg = 'Method not yet defined for mfd type %s' % mfd_type
        raise(mgs)

def combine_pt_sources(point_source_list, filename, name, nrml_version='04',
                       id_location_flag = 'location'):
    """Method for combining lists of point sources that are received
    and summing rates for co-located points
    Sources are joined based on id_location_flag, which can be 'id'
    or 'location' or None. Setting to None will mean at pts are 
    simply added together and no checking for co-located pts is undertake
    """
    # Get ids
    combined_pt_sources = []
    if id_location_flag is not None:
        for pt in point_source_list[0]:
            for source_model in point_source_list[1:]:
                #         print source_model
                for pt_source in source_model:
                    if id_location_flag == 'id':
                        if pt_source.source_id == pt.source_id:
                            new_rates = merge_rates(pt, pt_source)
                            pt.mfd.modify_set_mfd(pt.mfd.min_mag, pt.mfd.bin_width,
                                                  list(new_rates))
                            source_model.remove(pt_source)
                    elif id_location_flag == 'location':
                        #print type(pt_source)
                        #print type(pt)
                        # Check if location and nodal planes are the same
                        if pt_source.location.x == pt.location.x and \
                                pt_source.location.y == pt.location.y:
                            if pt_source.nodal_plane_distribution.data == pt.nodal_plane_distribution.data:
                                new_rates = merge_rates(pt, pt_source)
                                pt.mfd.modify_set_mfd(pt.mfd.min_mag, pt.mfd.bin_width,
                                                      list(new_rates))
                                source_model.remove(pt_source)
                        
    # once all overlapping point sources have been merged, add all to list
    # This should work as we have added rates to the first source model as
    # we have gone and removed sources in the same locations from the other
    # source mode lists
    for source_model in point_source_list:
        for pt in source_model:
            combined_pt_sources.append(pt)
    
    if nrml_version == '04':
#        for source in combined_pt_sources:
#            source_list.append(source)
#            id_index = max(id_index, source.source_id)
        nodes = list(map(obj_to_node, sorted(combined_pt_sources)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        with open(filename, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
    return combined_pt_sources
        

def write_combined_faults_points(point_sources, fault_sources,
                                 filename, name, area_sources = None,
                                 nrml_version='04'):
    """Write pts, area and fault sources to file
    :param point_sources:
       list without trt or dict with trt key of point sources
    """
    print 'Writing to source model file %s' % filename
    ps_id_index = 1
    fs_id_index = 1
    if nrml_version == '04':
        if type(point_sources) == dict:
            source_list = []
            for trt, sources in point_sources.iteritems():
                for source in sources:
                    source.source_id = 'PS_%i' % ps_id_index
                    source_list.append(source)
                    ps_id_index += 1
#                    id_index = max(id_index, source.source_id)
        elif type(point_sources) == list:
            source_list = copy.deepcopy(point_sources)
            for source in source_list:
                source.source_id = 'PS_%i' % ps_id_index
#                source_list.append(source)
                ps_id_index += 1
#                id_index = max(id_index, source.source_id)
        for fault_source in fault_sources:
#            id_index += 1
            fault_source.source_id = "FS_%i" % fs_id_index
            fs_id_index += 1
            source_list.append(fault_source)
        if area_sources is not None:
            for area_source in area_sources:
                source_list.append(area_source)
        nodes = list(map(obj_to_node, sorted(source_list)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        with open(filename, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
    elif nrml_version == '05':
        if type(point_sources) == dict:
            source_group_list = []
            id = 0
            for trt, sources in point_sources.iteritems():
                for source in sources:
                    id_index = max(id_index, source.source_id)
            for trt, sources in point_sources.iteritems():
                for fault_source in fault_sources:
                    if fault_source.tectonic_region_type == trt:
                        id_index += 1
                        fault_source.source_id = "%i" % id_index
                        sources.append(fault_source)
                if area_sources is not None:
                    for area_source in area_sources:
                        if area_source.tectonic_region_type == trt:
                            sources.append(area_source)
                source_group = SourceGroup(trt, sources = sources, id=id)
                id +=1
                source_group_list.append(source_group)
            write_source_model(filename, source_group_list,
                               name = name)
        elif type(point_sources) == list:
            msg = 'Method not yet implemented for nrml version 0.5'
            raise(msg)
    else:
        print 'Warning: nrml version not specfied, xml not created'

def read_simplefault_source(simplefault_source_file, rupture_mesh_spacing = 10):
    """Read nrml source model into simpmle fault objects
    """
    converter = SourceConverter(50, rupture_mesh_spacing, width_of_mfd_bin=0.1,
                                area_source_discretization=200.)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(simplefault_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(simplefault_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    for fault in sources:
       # print [method for method in dir(fault)]
       # print [method for method in dir(fault.mfd)]

#######Probably not actually needed now
        # Add max_mag attribute if needed:
#        try:
#            print fault.mfd.max_mag
#        except AttributeError:
#            min_mag, max_mag = fault.mfd.get_min_max_mag()
#            fault.mfd.max_mag = max_mag
#            print fault.mfd.max_mag
        pass
    return sources

def pt2fault_distance(pt_sources, fault_sources, min_distance = 5.0,
                      filename = 'source_model.xml',
                      buffer_distance = 100., nrml_version = '04',
                      name=None):

    """Calculate distances from a pt source rupture plane
    to the fault sources to then reduce Mmax on events that are 
    within a certain distance
    :param pt_sources:
        list of PointSource objects
    :param fault_sources:
        List of FaultSource objects
    :param min_distance:
        Minimum distance (km) within which we want a point source 
        rupture to be from a fault.
    :param filename:
        Name of output nrml file for revised pt source model
    :param buffer_distance:
        Km, initial filter to only process pts within this
        distance from the fault
    """

    if name is None:
        name = filename[:-4] + '_geom_filtered'
    id_index = 0 # We need to re-number all sources to avoid duplicate ids
    # Extract the points of the fault source mesh
    fault_lons = []
    fault_lats = []
    fault_depths = []
    for fault in fault_sources:
        whole_fault_surface = SimpleFaultSurface.from_fault_data(
            fault.fault_trace, fault.upper_seismogenic_depth,
            fault.lower_seismogenic_depth, fault.dip, 
            fault.rupture_mesh_spacing)
        fault_lons.append(whole_fault_surface.mesh.lons.flatten())
        fault_lats.append(whole_fault_surface.mesh.lats.flatten())
        fault_depths.append(whole_fault_surface.mesh.depths.flatten())
    fault_lons = np.concatenate(fault_lons)
    fault_lats = np.concatenate(fault_lats)
    fault_depths = np.concatenate(fault_depths)
    min_fault_lon = np.min(fault_lons)
    max_fault_lon = np.max(fault_lons)
    min_fault_lat = np.min(fault_lats)
    max_fault_lat = np.max(fault_lats)

    # Generate ruptures for point sources
    minimum_distance_list = []
    revised_point_sources = {'Cratonic': [], 'Non_cratonic': [], 
                             'Extended': [], 'Subduction': []}
    for pt in pt_sources:
        print 'Looping over point sources'
        # For speeding things up filter based on initial distances
        # to find points very far from or very close to a fault
        mfd_type = type(pt.mfd).__name__
        pt_depths = []
        for probs, depths in pt.hypocenter_distribution.data:
            pt_depths.append(depths)
        np_probs = []
        np_list = []
        for prob, nodal_plane in pt.nodal_plane_distribution.data:
            np_probs.append(prob)
            np_list.append(nodal_plane)
        centroid_distances = []
        for pt_depth in pt_depths:
            centroid_distances.append(distance(pt.location.longitude, pt.location.latitude,
                                      pt_depth, fault_lons, fault_lats, fault_depths))
        centroid_distances = np.array(centroid_distances).flatten()
  #      print 'Minimum distance', min(centroid_distances)
  #      print 'Maximum distance', max(centroid_distances)
        if (min(centroid_distances)) > buffer_distance:
            # Keep point as it, not within buffer distance of any faults
            revised_point_sources[pt.tectonic_region_type].append(pt)
            continue
        if (min(centroid_distances)) < min_distance:
            # Discard point sources as too close to a fault
            print 'Discarding point source, too close to a fault'
            continue
        rupture_mags = []
        rupture_lons = []
        rupture_lats = []
        rupture_depths = []
        rupture_strikes = []
        rupture_dips = []
        ruptures = pt.iter_ruptures()
        for rupture in ruptures:
            rupture_mags.append(rupture.mag)
            rupture_lons.append(rupture.surface.corner_lons)
            rupture_lats.append(rupture.surface.corner_lats)
            rupture_depths.append(rupture.surface.corner_depths)
            rupture_strikes.append(rupture.surface.strike)
            rupture_dips.append(rupture.surface.dip)
        rupture_mags = np.array(rupture_mags).flatten()
        # make the same length as the corners
        rupture_mags = np.repeat(rupture_mags, 4)
        rupture_strikes = np.repeat(rupture_strikes, 4)
        rupture_dips = np.repeat(rupture_dips, 4)
        rupture_lons = np.array(rupture_lons).flatten()
        rupture_lats = np.array(rupture_lats).flatten()
        rupture_depths = np.array(rupture_depths).flatten()
        print 'Doing meshgrid'
        lons1,lons2 = np.meshgrid(fault_lons, rupture_lons)
        lats1,lats2 = np.meshgrid(fault_lats, rupture_lats)
        depths1, depths2 = np.meshgrid(fault_depths, rupture_depths)

        # Calculate distance from pt to all fault
        print 'Distance calculations'
        distances = distance(lons1, lats1, depths1, lons2, lats2, depths2)
        closest_distance_to_faults = np.min(distances)
        print 'Shortest pt to fault distance is', closest_distance_to_faults
        minimum_distance_list.append(closest_distance_to_faults)

        # Find where the distance is less than the threshold min_distance
        too_close_lons = lons2[np.where(distances < min_distance)]
        too_close_lats = lats2[np.where(distances < min_distance)]
        if too_close_lons.size > 0:
            lon_indices = np.where(np.in1d(rupture_lons, too_close_lons))[0]
            lat_indices = np.where(np.in1d(rupture_lats, too_close_lats))[0]
            too_close_mags = rupture_mags[np.intersect1d(
                    lon_indices, lat_indices)]
            too_close_strikes = rupture_strikes[np.intersect1d(
                    lon_indices, lat_indices)]
            too_close_dips = rupture_dips[np.intersect1d(
                    lon_indices, lat_indices)]
        #    print 'Magnitudes of rupture close to fault', too_close_mags
        #    print 'Strikes of rupture close to fault', too_close_strikes
        #    print 'Dips of rupture close to fault', too_close_dips
            unique_strikes = np.unique(rupture_strikes)
            unique_dips = np.unique(rupture_dips)
            src_name_index = 0
            for prob, nodal_plane in pt.nodal_plane_distribution.data:
                id_index += 1
                src_name_index += 1
                # We are now splitting the source into many with different 
                # combinations of Mmaxs and nodal planes
                new_pt = copy.deepcopy(pt)
                new_pt.source_id ="%i" % id_index
                new_pt.name = new_pt.name + ("_%i" % src_name_index)
                new_np = NodalPlane(nodal_plane.strike, nodal_plane.dip, nodal_plane.rake)
                new_np_distribution = PMF([(1.0, new_np)]) # weight of nodal plane is 1 as making 
                # a separate source
                # Calculate new rates based on probability of original nodal plane
                new_pt.nodal_plane_distribution = new_np_distribution
                if mfd_type == 'TruncatedGRMFD':
                    b_val = pt.mfd.b_val
                    # rescale a value in log sapce
                    a_val = np.log10(np.power(10, pt.mfd.a_val)*prob)#*area_src_weight))
                    new_pt.mfd.modify_set_ab(a_val, b_val)
                elif mfd_type == 'EvenlyDiscretizedMFD':
                    mag_bins, rates = zip(*pt.mfd.get_annual_occurrence_rates())
                    mag_bins = np.array(mag_bins)
                    rates = np.array(rates)
                    new_rates = rates*prob#*area_src_weight)
                    new_pt.mfd.modify_set_mfd(new_pt.mfd.min_mag, new_pt.mfd.bin_width,
                                              list(new_rates))
                else:
                    msg = 'Weighting method for mfd type %s not yet defined' % mfd_type
                    raise(msg)
                pair_index = np.where(np.logical_and(too_close_strikes == nodal_plane.strike,
                                                         too_close_dips == nodal_plane.dip))
                 # Deal with intersecting cases
                if len(pair_index[0]) > 0:
                    intersecting_magnitudes = too_close_mags[pair_index]
                    minimum_magnitude_intersecting_fault = min(intersecting_magnitudes)
                    if minimum_magnitude_intersecting_fault >= \
                            (pt.mfd.min_mag + pt.mfd.bin_width):
                        new_mmax = minimum_magnitude_intersecting_fault - \
                                pt.mfd.bin_width
                        if mfd_type == 'TruncatedGRMFD':
                            new_pt.mfd.max_mag = new_mmax
                        if mfd_type == 'EvenlyDiscretizedMFD':
                            trimmed_rates = new_rates[np.where(mag_bins <= new_mmax)]
                    else:
                        print 'Minimum magnitude intersects fault, discarding source'
                        continue
                            
                else:
                    pass
                # Append revised source for given nodal plane distribution to 
                # list of revised sources
                print 'Appending revised source'
                revised_point_sources[pt.tectonic_region_type].append(new_pt)
        else:
            id_index += 1
            pt.source_id = "%i" % id_index
            'Appending original source'
            revised_point_sources[pt.tectonic_region_type].append(pt)
    if len(minimum_distance_list) > 0:
        print 'Overall minimum distance (km):', min(minimum_distance_list)

    # Write pts to source model on their own
    source_model_file = filename 
    print 'Writing to source model file %s' % source_model_file 
    if nrml_version == '04':
        source_list = []
        for trt, sources in revised_point_sources.iteritems():
            for source in sources:
                source_list.append(source)
        nodes = list(map(obj_to_node, sorted(source_list)))
        source_model = Node("sourceModel", {"name": name}, nodes=nodes)
        with open(source_model_file, 'wb') as f:
            nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
    elif nrml_version == '05':
        source_group_list = []
        id = 0
        for trt, sources in revised_point_sources.iteritems():
            source_group = SourceGroup(trt, sources = sources, id=id)
            id +=1
            source_group_list.append(source_group)
        write_source_model(source_model_file, source_group_list,
                           name = name)
    else:
        print 'Warning: nrml version not specfied, xml not created'

    # Write pts to source model with faults
    source_model_file = filename[:-4] +'_inc_faults.xml'
    name = name +'_inc_faults'
    write_combined_faults_points(revised_point_sources, fault_sources,
                                source_model_file, name, nrml_version='04')

def plot_sources(point_source, fault_source):
    from hmtk.plotting.mapping import HMTKBaseMap
    llon, ulon, llat, ulat = 105, 155, -45, -5,
    map_config = {'min_lon': np.floor(llon), 'max_lon': np.ceil(ulon),
                  'min_lat': np.floor(llat), 'max_lat': np.ceil(ulat), 'resolution':'i'}
    basemap1 = HMTKBaseMap(map_config, 'Point and fault sources')
    for pt_source in point_source:
         x,y = basemap1.m(pt_source.location.longitude, 
                         pt_source.location.latitude)
         basemap1.m.plot(x, y, 'rs',
                         markersize = 2.0)
    for simplefault in fault_source:
        trace_lons = np.array([pnt.longitude
                               for pnt in simplefault.fault_trace.points])
        trace_lats = np.array([pnt.latitude
                               for pnt in simplefault.fault_trace.points])
        x, y = basemap1.m(trace_lons, trace_lats)
        basemap1.m.plot(x, y, 'b', linewidth=1.3)
    basemap1.savemap('Point and fault sources.png')

if __name__ == "__main__":
    path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_inputs/NSHM'
    path = '/media/sf_openquake_shared_files/Australia/eq_hazard_tools/catalogue/'
    pt_source_file = os.path.join(path, 'source_model_leonard_2008_pts.xml')
    fault_path = '/home/openquake/GEM/eq_hazard_tools/shapefile2nrml/Adelaide_test/'
    simplefault_source_file = os.path.join(fault_path, 'Adelaide_faults.xml')
    pt_sources = read_pt_source(pt_source_file)
    fault_sources = read_simplefault_source(simplefault_source_file)
    plot_sources = plot_sources(pt_sources, fault_sources)
    revised_point_source_file = pt_source_file[:-4] + \
                                '_filter_near_fault_sources.xml'
    pt2fault_distance(pt_sources, fault_sources, 
                      filename = revised_point_source_file,
                      buffer_distance = 5.)
