"""Reads an input point source model and generates ruptures according 
to the MFD and nodal plane distirbution. Nearest distance to the fault is calculated and used to reduce the Mmax of the point source model. 
"""

import os, sys
import numpy as np
from openquake.commonlib.source import SourceModelParser
from openquake.commonlib.sourceconverter import SourceConverter
from openquake.commonlib.sourcewriter import write_source_model
from openquake.commonlib.node import Node
from openquake.commonlib.sourcewriter import obj_to_node
from openquake.commonlib import nrml
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
from openquake.hazardlib.geo.geodetic import distance, geodetic_distance
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
#plotting
from matplotlib import pyplot
from hmtk.plotting.mapping import HMTKBaseMap


def read_pt_source(pt_source_file):
    """Read nrml source model into pt source objects
    """
    converter = SourceConverter(50, 10, width_of_mfd_bin=0.1,
                                area_source_discretization=200.)
    parser = SourceModelParser(converter)
    try:
        sources = parser.parse_sources(pt_source_file)
    except AttributeError: # Handle version 2.1 and above
        sources = []
        groups = parser.parse_src_groups(pt_source_file)
        for group in groups:
            for source in group:
                sources.append(source)
    for pt in sources:
        print pt.mfd.max_mag
    return sources

def read_simplefault_source(simplefault_source_file):
    """Read nrml source model into simpmle fault objects
    """
    converter = SourceConverter(50, 2, width_of_mfd_bin=0.1,
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
        print fault.mfd.max_mag
    return sources

def pt2fault_distance(pt_sources, fault_sources, min_distance = 5):
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
    """

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
#    print fault_lons, type(fault_lons)
    # Generate ruptures for point sources
 
    minimum_distance_list = []
    buffer_distance = 4 # Degrees, initial filter to only process pts
    # within the region where we have faults
    for pt in pt_sources:
        # For speeding things up
        if pt.location.longitude < min_fault_lon - buffer_distance or \
           pt.location.longitude > max_fault_lon + buffer_distance or \
           pt.location.latitude < min_fault_lat - buffer_distance or \
           pt.location.latitude > max_fault_lat + buffer_distance:
            continue
#        print pt.location
        rupture_mags = []
        rupture_lons = []
        rupture_lats = []
        rupture_depths = []
        ruptures = pt.iter_ruptures()
        for rupture in ruptures:
            rupture_mags.append(rupture.mag)
            rupture_lons.append(rupture.surface.corner_lons)
            rupture_lats.append(rupture.surface.corner_lats)
            rupture_depths.append(rupture.surface.corner_depths)
        rupture_mags = np.array(rupture_mags).flatten()
        rupture_lons = np.array(rupture_lons).flatten()
        rupture_lats = np.array(rupture_lats).flatten()
        rupture_depths = np.array(rupture_depths).flatten()
        #Tile to get grid of distances
#        print np.shape(rupture_lons)
#        print np.shape(fault_lons)
        lons1,lons2 = np.meshgrid(fault_lons, rupture_lons)
        lats1,lats2 = np.meshgrid(fault_lats, rupture_lats)
        depths1, depths2 = np.meshgrid(fault_depths, rupture_depths)
#        print lats1[0]
#        print lats2[0]
         # Calculate distance from pt to all faults, then
        # find the minimum and equate to magnitude based on
        # scaling relationship
        distances = distance(lons1, lats1, depths1, lons2, lats2, depths2)
 #       print 'distances', distances
        closest_distance_to_faults = np.min(distances)
        print closest_distance_to_faults
        minimum_distance_list.append(closest_distance_to_faults)
    print 'Overall minimum', min(minimum_distance_list)
#        pyplot.clf()
#        pyplot.scatter(lons1, lons2)
#        pyplot.savefig('lons.png')
#        sys.exit()
        

def plot_sources(point_source, fault_source):
   
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
    pt2fault_distance(pt_sources, fault_sources)
