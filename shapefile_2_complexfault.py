import os
import ogr, gdal

def parse_line_shapefile(shapefile, shapefile_depth_attribute):
    """Read the line shapefile with contours on the fault plain

        Return a list of lists of [x,y,z] points defining each line
    """


    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile,0)

    layer = data_source.GetLayer()

    fault_contours = []

    for feature in layer:
        line = feature.GetGeometryRef().GetPoints()
        depth = float(feature.GetField(shapefile_depth_attribute))
        line = [list(pts)+[depth] for pts in line]
        fault_contours.append(line)

    msg = 'Line shapefile must contain at least 2 fault source contour lines'
    assert len(fault_contours)>1, msg

    # Make sure the ordering of depth is from smallest to largest
    new_fault_contours = []
    contour_depths = [ fc[0][2] for fc in fault_contours ]
    contour_depths.sort()

    for cd in contour_depths:
        for j in range(len(fault_contours)):
            if fault_contours[j][0][2] == cd:
                new_fault_contours.append(fault_contours[j])
                break

    return new_fault_contours



def append_xml_header(output_xml, 
                      source_model_name, 
                      complex_fault_id, 
                      complex_fault_name, 
                      complex_fault_tectonic_region):
    """Add information to the nrml which goes above the fault contours

    """

    # First line = blank
    output_xml.append('')

    # Various header info which the user does not control
    output_xml.append("<?xml version='1.0' encoding='utf-8'?>")
    output_xml.append('<nrml xmlns:gml="http://www.opengis.net/gml"' + \
        ' xmlns="http://openquake.org/xmlns/nrml/0.4">')
   
    # Source model information 
    output_xml.append('  <sourceModel name=' + str(source_model_name) + '>')
    output_xml.append('  <complexFaultSource id="' + str(complex_fault_id) + \
        '"' +  ' name="' + str(complex_fault_name) + '"' + \
        ' tectonicRegion="' + str(complex_fault_tectonic_region) +'">')

    output_xml.append('')

    return


def append_gml_Linestring(output_xml, fc):
    """Convenience function to append the xyz coordinates as a gml Linestring
        
        @param output_xml List holding the lines of the output xml (which is being built)
        @param fc List of lists of the form [x,y,z], defining the contour

        @return nothing, but output_xml is modified with the gml Linestring appended 
    
    """
    output_xml.append('          <gml:Linestring>')
    output_xml.append('            <gml:posList>')

    # Add the geometry 
    for i in range(len(fc)):
        output_xml.append('               ' + str(fc[i][0]) + ' ' + str(fc[i][1]) + ' ' + str(fc[i][2]) )

    # Footer
    output_xml.append('            <gml:posList>')
    output_xml.append('          <gml:Linestring>')

    return


def append_rupture_geometry(output_xml, fault_contours):
    """Append the fault contours to the nrml

    """

    ## Top edge     

    # Header
    output_xml.append('      <complexFaultGeometry>')
    output_xml.append('        <faultTopEdge>')
    append_gml_Linestring(output_xml, fault_contours[0])
    output_xml.append('        <faultTopEdge>')
    output_xml.append('')

    # Intermediate edges
    if len(fault_contours)>2:
        for i in range(1, len(fault_contours)-1):
            output_xml.append('        <intermediateEdge>')
            append_gml_Linestring(output_xml, fault_contours[i])
            output_xml.append('        <intermediateEdge>')
            output_xml.append('')
           
    
    # Bottom edge 
    output_xml.append('        <faultBottomEdge>')
    append_gml_Linestring(output_xml, fault_contours[-1])
    output_xml.append('        <faultBottomEdge>')
    output_xml.append('      <complexFaultGeometry>')
    output_xml.append('')

    return


def append_earthquake_information(output_xml, magnitude_scaling_relation,
        rupture_aspect_ratio, a_value, b_value, min_mag, max_mag, rake,):
    """

    """
    output_xml.append('      <magScaleRel><' + \
        str(magnitude_scaling_relation) + '>')
    output_xml.append('')
      
    output_xml.append('      <ruptAspectRatio><' + str(rupture_aspect_ratio) + '>')
    output_xml.append('')

    output_xml.append('      <truncGutenbergRichterMFD aValue=' + \
        str(a_value) + ' bValue=' + str(b_value) +' minMag=' + str(min_mag) + \
        ' maxMag=' + str(max_mag) + '/>')
    output_xml.append('')

    output_xml.append('      <rake><' + str(rake) + '>')
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

    append_xml_header(output_xml, source_model_name, complex_fault_id, 
        complex_fault_name, complex_fault_tectonic_region)
    
    append_rupture_geometry(output_xml, fault_contours)

    append_earthquake_information(output_xml, magnitude_scaling_relation,
        rupture_aspect_ratio, a_value, b_value, min_mag, max_mag, rake) 

    # Add newlines
    output_xml = [ oxml+'\n' for oxml in output_xml ]

    return output_xml

############################################################################

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Convert a line shapefile to an openquake complex fault in nrml format')

    parser.add_argument('-shapefile', type=str, default="",
        help='Filename of line shapefile with rupture plain contours. Each line geometry has an attribute defining its depth')

    parser.add_argument('-shapefile_depth_attribute', type=str, default=None,
        help='Name of the attribute table column (in shapefile) which contains the depth of each contour')

    parser.add_argument('-source_model_name', type=str,
        default='unnamed_source', help='Name for the source model')

    parser.add_argument('-complex_fault_id', type=int, default=0,
        help='integer id for the complex_fault')
    
    parser.add_argument('-complex_fault_name', type=str, 
        default='unnamed_complex_fault', help='Name for the complex fault')
    
    parser.add_argument('-complex_fault_tectonic_region', type=str, 
        default='unnamed_complex_fault_region',
        help='Name for the tectonic region of the complex fault')

    parser.add_argument('-magnitude_scaling_relation', type=str, 
        default='/magScaleRel',
        help='Openquake magnitude scaling relation descriptor (see manual)')
    
    parser.add_argument('-rupture_aspect_ratio', type=str, 
        default='/ruptAspectRatio',
        help='Openquake rupture aspect ratio (see manual)')
    
    parser.add_argument('-a_value', type=str, 
        default='""',
        help='Openquake truncated Gutenberg Richter MFD a-value (see manual)')
    
    parser.add_argument('-b_value', type=str, 
        default='""',
        help='Openquake truncated Gutenberg Richter MFD b-value (see manual)')
    
    parser.add_argument('-min_mag', type=str, 
        default='""',
        help='Openquake truncated Gutenberg Richter MFD minMag (see manual)')
    
    parser.add_argument('-max_mag', type=str, 
        default='""',
        help='Openquake truncated Gutenberg Richter MFD maxMag (see manual)')
    
    parser.add_argument('-rake', type=str, 
        default='/rake',
        help='earthquake rake')

    parser.add_argument('-output_dir', type=str, default='nrml',
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

    f = open( os.path.join(args.output_dir, args.source_model_name + '.xml'), 
            'w')
    f.writelines(output_xml_text)
    f.close()

