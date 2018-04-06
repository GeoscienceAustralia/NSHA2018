"""

Command line tool to convert line shapefile + various input parameters to an
openquake complex fault in nrml format

Run with 
> python shapefile_2_complexfault.py -h
to get help

Gareth Davies, Geoscience Australia, 2014

"""

import os
import ogr
import argparse
from numpy import sqrt
import numpy
from openquake.hazardlib.geo.utils import spherical_to_cartesian

def check_aki_richards_convention(edges):
    # Adapted from openquake/hazardlib/geo/surface/complex_fault.py
    ul, ur, bl, br = spherical_to_cartesian(
        [edges[0][0][0], edges[0][-1][0], edges[-1][0][0], edges[-1][-1][0]],
        [edges[0][0][1], edges[0][-1][1], edges[-1][0][1], edges[-1][-1][1]],
        [edges[0][0][2], edges[0][-1][2], edges[-1][0][2], edges[-1][-1][2]])
    top_edge = ur - ul
    left_edge = bl - ul
    right_edge = br - ur
    left_cross_top = numpy.cross(left_edge, top_edge)
    right_cross_top = numpy.cross(right_edge, top_edge)
    
    left_cross_top /= numpy.sqrt(numpy.dot(left_cross_top, left_cross_top))
    right_cross_top /= numpy.sqrt(
        numpy.dot(right_cross_top, right_cross_top)
        )
    ul /= numpy.sqrt(numpy.dot(ul, ul))
    ur /= numpy.sqrt(numpy.dot(ur, ur))
    
    # rounding to 1st digit, to avoid ValueError raised for floating point                                                                                    
    # imprecision                                                                                                                                             
    angle_ul = round(
        numpy.degrees(numpy.arccos(numpy.dot(ul, left_cross_top))), 1
        )
    angle_ur = round(
        numpy.degrees(numpy.arccos(numpy.dot(ur, right_cross_top))), 1
        )

    if (angle_ul > 90) or (angle_ur > 90):
        print 'Contours do not conform to Aki-Richards convention, re-ordering'
        new_contours = []
        for contour in edges:
            new_contours.append(contour[::-1])
        edges=new_contours
    return edges

def parse_line_shapefile(shapefile, shapefile_depth_attribute,
                         boundary=[-360, 360, -90, 90]):
    """Read the line shapefile with contours on the fault plain
    Bounday: [min_lon, max_lon, min_lat, max_lat] Boundary to clip fault
    segment by 
        Return a list of lists of [x,y,z] points defining each line
    """

    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile, 0)

    layer = data_source.GetLayer()

    fault_contours = []

    m=0
    for feature in layer:
        line = feature.GetGeometryRef().GetPoints()
        depth = float(feature.GetField(shapefile_depth_attribute))
        pt_list = [list(pts) + [depth] for pts in line]
        lons = []
        lats = []
        depths = []
        for pt in pt_list:
            lons.append(pt[0])
            lats.append(pt[1])
            depths.append(pt[2])
        if m==0:
            start_lon=lons[0]
            start_lat=lats[0]
        if m > 0:
            # Check if start point is closer to start of top trace than
            # end point; if not, re-order to ensure all contours are
            # in the same order. Note we assume a flat earth, 
            # which should be ok in most circumstances as errors should be
            # 'small' relative to fault lengths
            start_dist = sqrt((start_lon-lons[0])**2 + (start_lat-lats[0])**2)
            end_dist = sqrt((start_lon-lons[-1])**2 + (start_lat-lats[-1])**2)
            if end_dist < start_dist:
                print 'reordering points to make contours go in same direction'
                lons = lons[::-1]
                lats = lats[::-1]
        # Check Aki-Richards convention
        
        if boundary is not None:
            line = [ (i,j,k) for (i,j,k) in zip(lons,lats,depths) if i >= boundary[0] if i <= boundary[1] if j >= boundary[2] if j <= boundary[3]]
        else:
            line = [(i,j,k) for (i,j,k) in zip(lons,lats,depths)]
        m+=1
        fault_contours.append(line)
    # Check Aki-Richards convention 
    new_contours = check_aki_richards_convention(fault_contours)
    # And double-check
    if new_contours != fault_contours:
        new_contours = check_aki_richards_convention(new_contours)
    fault_contours=new_contours
    msg = 'Line shapefile must contain at least 2 fault source contour lines'
    assert len(fault_contours) > 1, msg

    # Make sure the ordering of depth is from smallest to largest
    new_fault_contours = []
    contour_depths = [fc[0][2] for fc in fault_contours]
    contour_depths.sort()

    for cd in contour_depths:
        for j in range(len(fault_contours)):
            if fault_contours[j][0][2] == cd:
                new_fault_contours.append(fault_contours[j])
                break

    return new_fault_contours


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
    return

def append_fault_source_header(output_xml,
                               complex_fault_id,
                               complex_fault_name,
                               complex_fault_tectonic_region):
    output_xml.append('  <complexFaultSource id="' + str(complex_fault_id) +
                      '"' + ' name="' + str(complex_fault_name) + '"' +
                      ' tectonicRegion="' +
                      str(complex_fault_tectonic_region) + '">')

    output_xml.append('')

    return


def append_gml_Linestring(output_xml, fc, dps=10):
    """Convenience function to append the xyz coordinates as a gml Linestring

        @param output_xml List holding lines of the output xml being built
        @param fc List of lists of the form [x,y,z], defining the contour
        @param nps Number of decimal points to limit output coordinates to
        @return nothing, but the gml linestring is appended to the output_xml

    """
    output_xml.append('          <gml:LineString>')
    output_xml.append('            <gml:posList>')

    # Add the geometry
    for i in range(len(fc)):
        try:
            output_xml.append(
                '               ' + '{:.{prec}f}'.format(fc[i][0], prec=dps) + ' ' +\
                    '{:.{prec}f}'.format(fc[i][1], prec=dps) + ' ' +\
                    str(fc[i][2]))
        except ValueError:
            print 'Needs python 2.7 or higher for string formatting'
            sys.exit()
    # Footer
    output_xml.append('            </gml:posList>')
    output_xml.append('          </gml:LineString>')

    return


def append_rupture_geometry(output_xml, fault_contours, dps=10):
    """Append the fault contours to the nrml
    dps = number of decimal places for coordinates 
    """

    # Top edge

    # Header
    output_xml.append('      <complexFaultGeometry>')
    output_xml.append('        <faultTopEdge>')
    append_gml_Linestring(output_xml, fault_contours[0], dps)
    output_xml.append('        </faultTopEdge>')
    output_xml.append('')

    # Intermediate edges
    if len(fault_contours) > 2:
        for i in range(1, len(fault_contours) - 1):
            output_xml.append('        <intermediateEdge>')
            append_gml_Linestring(output_xml, fault_contours[i], dps)
            output_xml.append('        </intermediateEdge>')
            output_xml.append('')

    # Bottom edge
    output_xml.append('        <faultBottomEdge>')
    append_gml_Linestring(output_xml, fault_contours[-1], dps)
    output_xml.append('        </faultBottomEdge>')
    output_xml.append('      </complexFaultGeometry>')
    output_xml.append('')

    return


def append_earthquake_information(output_xml, magnitude_scaling_relation,
                                  rupture_aspect_ratio, a_value, b_value,
                                  min_mag, max_mag, rake,):
    """

    """
    output_xml.append('      <magScaleRel>' +
                      str(magnitude_scaling_relation) + '</magScaleRel>')
    output_xml.append('')

    output_xml.append(
        '      <ruptAspectRatio>' + str(rupture_aspect_ratio) + '</ruptAspectRatio>')
    output_xml.append('')

    output_xml.append('      <truncGutenbergRichterMFD aValue="' +
                      str(a_value) + '" bValue="' + str(b_value) +
                      '" minMag="' + str(min_mag) +
                      '" maxMag="' + str(max_mag) + '" />')
    output_xml.append('')

    output_xml.append('      <rake>' + str(rake) + '</rake>')
    output_xml.append('    </complexFaultSource>')
    output_xml.append('  </sourceModel>')

    output_xml.append('</nrml>')

    return


def nrml_from_shapefile(shapefile,
                        shapefile_depth_attribute,
                        source_model_name,
                        complex_fault_id,
                        complex_fault_name,
                        complex_fault_tectonic_region,
                        magnitude_scaling_relation,
                        rupture_aspect_ratio,
                        a_value,
                        b_value,
                        min_mag,
                        max_mag,
                        rake,
                        output_dir,
                        quiet):
    """Driver routine to convert nrml to shapefile

    """
    # Get geometry
    fault_contours = parse_line_shapefile(shapefile, shapefile_depth_attribute)

    # Output is written line-by-line to this list
    output_xml = []

    append_xml_header(output_xml, source_model_name)

    append_fault_source_header(output_xml, complex_fault_id,
                      complex_fault_name, complex_fault_tectonic_region)

    append_rupture_geometry(output_xml, fault_contours)

    append_earthquake_information(output_xml, magnitude_scaling_relation,
                                  rupture_aspect_ratio, a_value, b_value,
                                  min_mag, max_mag, rake)

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

    parser.add_argument('-shapefile_depth_attribute', type=str, default=None,
                        help='Name of the attribute table column ' +
                        '(in shapefile) which contains the depth of each ' +
                        'contour')

    parser.add_argument('-source_model_name', type=str,
                        default='unnamed_source',
                        help='Name for the source model')

    parser.add_argument('-complex_fault_id', type=int, default=0,
                        help='integer id for the complex_fault')

    parser.add_argument('-complex_fault_name', type=str,
                        default='unnamed_complex_fault',
                        help='Name for the complex fault')

    parser.add_argument(
        '-complex_fault_tectonic_region',
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
                                          args.shapefile_depth_attribute,
                                          args.source_model_name,
                                          args.complex_fault_id,
                                          args.complex_fault_name,
                                          args.complex_fault_tectonic_region,
                                          args.magnitude_scaling_relation,
                                          args.rupture_aspect_ratio,
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
