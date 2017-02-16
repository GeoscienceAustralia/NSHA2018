"""

Command line tool to convert line shapefile + various input parameters to an
openquake simple fault in nrml format

Shapefile may contain multiple lines respresenting multiple 
fault traces

Run with 
> python shapefile_2_simplefault.py -h
to get help

Gareth Davies, Geoscience Australia, 2014
Jonathan Griffin, Geoscience Australia, June 2016

"""

import os
import argparse
import ogr
import shapefile
from geopy import distance
from shapely.geometry import Point, Polygon
import numpy as np
from NSHA2018.mfd import fault_slip_rate_GR_conversion
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR

def parse_line_shapefile(shapefile,shapefile_faultname_attribute,
                         shapefile_dip_attribute, 
                         shapefile_sliprate_attribute,
                         shapefile_uplift_attribute):
    """Read the line shapefile with of fault surface traces and
    extrace data from attibute table
        Return a list of lists of [x,y] points defining each line
    and calculates the length of the fault
    """

    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile, 0)

    layer = data_source.GetLayer()

    fault_traces = []
    faultnames = []
    dips = []
    sliprates = []
    fault_lengths = []
    distance.VincentyDistance.ELLIPSOID = 'WGS_84'
    d = distance.distance
    for feature in layer:
        line = feature.GetGeometryRef().GetPoints()
        try:
            faultname = str(feature.GetField(shapefile_faultname_attribute))
        except ValueError:
            faultname = '""'
        faultnames.append(faultname)
        try:
            dip = float(feature.GetField(shapefile_dip_attribute))
            # Assume dips=0 is wrong, set to 45 degreee
            if dip == 0.0:
                dip = 45.0
        except ValueError:
            dip = '""'
        dips.append(dip)
        try:
            sliprate = float(feature.GetField(shapefile_sliprate_attribute))
            # Convert from m/ma to mm/a
            sliprate = sliprate/1000
        except ValueError:
            sliprate = '""'
        # If sliprate not given, calculate from uplift rate
        if sliprate == "" or sliprate == 0.0:
            try:
                upliftrate = float(feature.GetField(shapefile_uplift_attribute))
                # Convert from m/ma to mm/a
                upliftrate = upliftrate/1000
                # Calculate sliprate using dip and uplift rate
                sliprate = upliftrate*np.tan(dip*np.pi/180.)
            except ValueError:
                sliprate = '""'
        sliprates.append(sliprate)
        line = [list(pts) for pts in line]
        fault_length = 0
        for i in range(len(line)):
            if i == len(line) -1:
                break
            else:
                pt_distance = d(line[i],line[i+1]).meters/1000
                fault_length+=pt_distance
        fault_lengths.append(fault_length)
        fault_traces.append(line)

    msg = 'Line shapefile must contain at least 1 fault'
    assert len(fault_traces) > 0, msg
#    print fault_traces
#    print len(fault_traces)

    # Make sure the ordering of depth is from smallest to largest
#    new_fault_contours = []
#    contour_depths = [fc[0][2] for fc in fault_contours]
#    contour_depths.sort()

#    for cd in contour_depths:
#        for j in range(len(fault_contours)):
#            if fault_contours[j][0][2] == cd:
#                new_fault_contours.append(fault_contours[j])
#                break


    return fault_traces, faultnames, dips, sliprates, fault_lengths

def b_value_from_region(fault_traces, region_shapefile):
    """Get regional b-values for each fault
    """
    print 'Getting b-value from Leonard2008...'

    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(region_shapefile, 0)
    lsf = data_source.GetLayer()
    lbval = []
    for feature in lsf:
        lbval.append(float(feature.GetField('BVAL_BEST')))
#    geom = feature.GetGeometryRef()
    lsf = shapefile.Reader(region_shapefile)
    # get Leonard polygons                                      
    l08_shapes = lsf.shapes()
    # loop through faults and find point in poly
    b_values = []
    for fault_trace in fault_traces:
        trace_b_list = []
        for zone_bval, l_shape in zip(lbval, l08_shapes):
            l_poly = Polygon(l_shape.points)
        # check if leonard centroid in domains poly
            for point in fault_trace:
                pt = Point(point[0], point[1])
                if pt.within(l_poly):
                    bval = float(zone_bval)
                    trace_b_list.append(bval)
        # Find most common b_value
        try:
            (values,counts) = np.unique(trace_b_list,return_counts=True)
            ind=np.argmax(counts)
            b_values.append(values[ind])
        except ValueError:
            b_values.append(0.85) # Default value for points outside Leonard zones
   # print b_values
    return b_values

def trt_from_domains(fault_traces, domains_shapefile):
    """Get tectonic region type from domains
    """

    print 'Getting tectonic region type from Domains shapefile'
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(domains_shapefile, 0)
    dsf = data_source.GetLayer()
    trt_types = []
    for feature in dsf:
        trt_types.append(feature.GetField('TRT'))
#    geom = feature.GetGeometryRef()
    dsf = shapefile.Reader(domains_shapefile)
    # get Leonard polygons                                                                  
    dom_shapes = dsf.shapes()
    # loop through faults and find point in poly                                            
    trt_list = []
    for fault_trace in fault_traces:
        trace_trt_list = []
        for zone_trt, dom_shape in zip(trt_types, dom_shapes):
            dom_poly = Polygon(dom_shape.points)
        # check if leonard centroid in domains poly                                         
            for point in fault_trace:
                pt = Point(point[0], point[1])
                if pt.within(dom_poly):
                    trt = zone_trt
                    trace_trt_list.append(trt)
        # Find most common trt
        try:
            (values,counts) = np.unique(trace_trt_list,return_counts=True)
            ind=np.argmax(counts)
            trt_list.append(values[ind])
        except ValueError:
            print 'Warning: setting fault TRT to default value'
            trt_list.append('Non_cratonic') # Default value for points outside Domains model
   # print b_values                                                
    return trt_list

def append_xml_header(output_xml,
                      source_model_name):
    """Add information to the nrml which goes above the fault contours

    """

    # I think this breaks openquake
    # First line = blank
#    output_xml.append('')

    # Various header info which the user does not control
    output_xml.append("<?xml version='1.0' encoding='utf-8'?>")
    output_xml.append('<nrml xmlns:gml="http://www.opengis.net/gml"' +
                      ' xmlns="http://openquake.org/xmlns/nrml/0.4">')

    # Source model information
    output_xml.append('  <sourceModel name="' + str(source_model_name) + '">')
    #output_xml.append('')

    return


def append_gml_Linestring(output_xml, fc):
    """Convenience function to append the xyz coordinates as a gml Linestring

        @param output_xml List holding lines of the output xml being built
        @param fc List of lists of the form [x,y,z], defining the contour

        @return nothing, but the gml linestring is appended to the output_xml

    """
    output_xml.append('          <gml:LineString>')
    output_xml.append('            <gml:posList>')

    # Add the geometry
    for i in range(len(fc)):
        output_xml.append(
            '               ' + str(fc[i][0]) + ' ' + str(fc[i][1]))

    # Footer
    output_xml.append('            </gml:posList>')
    output_xml.append('          </gml:LineString>')

    return


def append_rupture_geometry(output_xml, trace, dip,
                            simple_fault_id, simple_fault_name,
                            upper_depth, lower_depth, 
                            simple_fault_tectonic_region):
    """Append the fault contours to the nrml

    """

    # Basic fault info
    output_xml.append('')
    output_xml.append('  <simpleFaultSource id="' + str(simple_fault_id) +
                      '"' + ' name="' + str(simple_fault_name) + '"' +
                      ' tectonicRegion="' +
                      str(simple_fault_tectonic_region) + '">')
    # Geometry
    output_xml.append('      <simpleFaultGeometry>')
    append_gml_Linestring(output_xml, trace)
#    output_xml.append('')
    output_xml.append('          <dip>' + str(dip) + '</dip>')
    output_xml.append('          <upperSeismoDepth>' + str(upper_depth) + '</upperSeismoDepth>')
    output_xml.append('          <lowerSeismoDepth>' + str(lower_depth) + '</lowerSeismoDepth>')
    output_xml.append('      </simpleFaultGeometry>')
#    output_xml.append('')
    return


def append_earthquake_information(output_xml, magnitude_scaling_relation,
                                  rupture_aspect_ratio, a_value, b_value,
                                  min_mag, max_mag, rake,):
    """

    """
    output_xml.append('      <magScaleRel>' +
                      str(magnitude_scaling_relation) + '</magScaleRel>')
  #  output_xml.append('')

    output_xml.append(
        '      <ruptAspectRatio>' + str(rupture_aspect_ratio) + '</ruptAspectRatio>')
 #   output_xml.append('')

    output_xml.append('      <truncGutenbergRichterMFD aValue="' +
                      str(a_value) + '" bValue="' + str(b_value) +
                      '" minMag="' + str(min_mag) +
                      '" maxMag="' + str(max_mag) + '" />')
#    output_xml.append('')

    output_xml.append('      <rake>' + str(rake) + '</rake>')
    output_xml.append('    </simpleFaultSource>')

    return


def nrml_from_shapefile(shapefile,
                        shapefile_faultname_attribute,
                        shapefile_dip_attribute,
                        shapefile_sliprate_attribute,
                        source_model_name,
                        simple_fault_tectonic_region,
                        magnitude_scaling_relation,
                        rupture_aspect_ratio,
                        upper_depth,
                        lower_depth,
                        a_value,
                        b_value,
                        min_mag,
                        max_mag,
                        rake,
                        output_dir,
                        shapefile_uplift_attribute=None,
                        quiet = True):
    """Driver routine to convert nrml to shapefile
     
    """
    # Get geometry
    fault_traces, faultnames, dips, \
    sliprate, fault_lengths = parse_line_shapefile(shapefile,
                                                   shapefile_faultname_attribute,
                                                   shapefile_dip_attribute, 
                                                   shapefile_sliprate_attribute,
                                                   shapefile_uplift_attribute)

     # Output is written line-by-line to this list
    output_xml = []

    append_xml_header(output_xml, source_model_name)

    # If b-value is not given, take from Leonard 2008 model
    region_shapefile = '../zones/Leonard2008/shapefiles/LEONARD08_NSHA18_MFD.shp'
    if b_value is None:
        b_value = b_value_from_region(fault_traces, region_shapefile)
    
    # If tectonic region type is not given, take from domains model
    domains_shapefile = '../zones/Domains/shapefiles/DOMAINS_NSHA18.shp'
    if simple_fault_tectonic_region is None:
        simple_fault_tectonic_region = trt_from_domains(fault_traces, domains_shapefile)

    # Loop through each fault and add source specific info
    for i in range(len(fault_traces)):
        # Skip faults with zero or null sliprate
        if sliprate[i] == "" or sliprate[i] == 0:
            continue
        simple_fault_id = i
        A = fault_lengths[i]*(float(lower_depth)-float(upper_depth))
        # Calculate M_max from scaling relations
        scalrel = Leonard2014_SCR()
        max_mag = scalrel.get_median_mag(A, float(rake))
#        print A
        # Calculate GR a values from slip rate
        if sliprate[i] != '""':
           # print sliprate[i]
            a_value, moment_rate = fault_slip_rate_GR_conversion.slip2GR(sliprate[i], A,
                                                                         float(b_value[i]), 
                                                                         float(max_mag),
                                                                         M_min=0.0)
        append_rupture_geometry(output_xml, fault_traces[i],
                                dips[i], simple_fault_id,
                                faultnames[i], upper_depth,
                                lower_depth, simple_fault_tectonic_region[i])

        append_earthquake_information(output_xml,
                                      magnitude_scaling_relation,
                                      rupture_aspect_ratio, a_value, b_value[i],
                                      min_mag, max_mag, rake)

    # Close xml
    output_xml.append('  </sourceModel>')
    output_xml.append('</nrml>')

    # Add newlines
    output_xml = [oxml + '\n' for oxml in output_xml]

    return output_xml

############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert a line shapefile to an openquake complex fault' +
        ' in nrml format')

    parser.add_argument('-shapefile', type=str, default="",
                        help='Filename of line shapefile with rupture ' +
                        'plain contours. Each line geometry has an ' +
                        'attribute defining its depth')

    parser.add_argument('-shapefile_faultname_attribute', type=str, default=None,
                        help='Name of the attribute table column ' +
                        '(in shapefile) which contains the name of each ' +
                        'fault')

    parser.add_argument('-shapefile_dip_attribute', type=str, default=None,
                        help='Name of the attribute table column ' +
                        '(in shapefile) which contains the dip of each ' +
                        'fault')

    parser.add_argument('-shapefile_sliprate_attribute', type=str, default=None,
                        help='Name of the attribute table column ' +
                        '(in shapefile) which contains the sliprate of each ' +
                        'fault')

    parser.add_argument('-source_model_name', type=str,
                        default='unnamed_source',
                        help='Name for the source model')

    parser.add_argument(
        '-simple_fault_tectonic_region',
        type=str,
        default='unnamed_complex_fault_region',
        help='Name for the tectonic region of the complex fault')

    parser.add_argument(
        '-magnitude_scaling_relation',
        type=str,
        default='/magScaleRel',
        help='Openquake magnitude scaling relation descriptor (see manual)')

    parser.add_argument('-rupture_aspect_ratio', type=str,
                        default='/ruptAspectRatio',
                        help='Openquake rupture aspect ratio (see manual)')
  
    parser.add_argument('-upper_depth', type=str,
                        default='""',
                        help='Openquake upper seismogenic depth (see manual)')

    parser.add_argument('-lower_depth', type=str,
                        default='""',
                        help='Openquake lower seismogenic depth (see manual)')

    parser.add_argument(
        '-a_value',
        type=str,
        default='""',
        help='Openquake truncated Gutenberg Richter MFD a-value (see manual)')

    parser.add_argument(
        '-b_value',
        type=str,
        default='""',
        help='Openquake truncated Gutenberg Richter MFD b-value (see manual)')

    parser.add_argument(
        '-min_mag',
        type=str,
        default='""',
        help='Openquake truncated Gutenberg Richter MFD minMag (see manual)')

    parser.add_argument(
        '-max_mag',
        type=str,
        default='""',
        help='Openquake truncated Gutenberg Richter MFD maxMag (see manual)')

    parser.add_argument('-rake', type=str,
                        default='/rake',
                        help='earthquake rake')

    parser.add_argument(
        '-output_dir',
        type=str,
        default='nrml',
        help='Location for output files. Created if necessary.')

    parser.add_argument('-quiet', action='store_true', default=False,
                        help="Don't print anything")

    args = parser.parse_args()

    try:
        shapefile_exists = os.path.exists(args.shapefile)
    except:
        parser.print_help()
        raise Exception('shapefile ' + str(args.shapefile) + ' not found')

    if not shapefile_exists:
        parser.print_help()
        raise Exception('shapefile ' + str(args.shapefile) + ' not found')

    # Main code here
    output_xml_text = nrml_from_shapefile(args.shapefile,
                                          args.shapefile_faultname_attribute,
                                          args.shapefile_dip_attribute,
                                          args.shapefile_sliprate_attribute,
                                          args.source_model_name,
                                          args.simple_fault_tectonic_region,
                                          args.magnitude_scaling_relation,
                                          args.rupture_aspect_ratio,
                                          args.upper_depth,
                                          args.lower_depth,
                                          args.a_value,
                                          args.b_value,
                                          args.min_mag,
                                          args.max_mag,
                                          args.rake,
                                          args.output_dir,
                                          args.quiet)

    # Write to file
    try:
        os.mkdir(args.output_dir)
    except:
        pass

    f = open(os.path.join(args.output_dir, args.source_model_name + '.xml'),
             'w')
    f.writelines(output_xml_text)
    f.close()
