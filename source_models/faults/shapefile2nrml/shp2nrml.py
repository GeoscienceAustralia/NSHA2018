"""Putting generic tools for extracting data from 
shapefiles into one place
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
from openquake.hazardlib.mfd import YoungsCoppersmith1985MFD, TruncatedGRMFD

def parse_line_shapefile(shapefile,shapefile_faultname_attribute,
                         shapefile_dip_attribute, 
                         shapefile_sliprate_attribute,
                         shapefile_uplift_attribute=None,
                         slip_units = 'm/ma'):
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
            if slip_units == 'm/ma':
                sliprate = sliprate/1000
            elif slip_units == 'mm/a':
                sliprate = sliprate
            else:
                raise Exception('Unkown sliprate units of %s' % slip_units)
        except ValueError:
            sliprate = '""'
        # If sliprate not given, calculate from uplift rate
        if sliprate == "" or sliprate == 0.0:
            if shapefile_uplift_attribute is not None:
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

    return fault_traces, faultnames, dips, sliprates, fault_lengths

def b_value_from_region(fault_traces, region_shapefile, default_b = 1.0):
    """Get regional b-values for each fault
    """
    print 'Getting b-value from zones in %s' % region_shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(region_shapefile, 0)
    lsf = data_source.GetLayer()
    lbval = []
    for feature in lsf:
        lbval.append(float(feature.GetField('BVAL_BEST')))
    lsf = shapefile.Reader(region_shapefile)
    l08_shapes = lsf.shapes()
    b_values = []
    for fault_trace in fault_traces:
        trace_b_list = []
        for zone_bval, l_shape in zip(lbval, l08_shapes):
            l_poly = Polygon(l_shape.points)
        # check if region centroid in domains poly
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
            b_values.append(default_b) # Default value for undefined points
    return b_values

def trt_from_domains(fault_traces, domains_shapefile, 
                     default_trt = 'Non_cratonic'):
    """Get tectonic region type from domains
    """
    print 'Getting tectonic region type from %s' % domains_shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(domains_shapefile, 0)
    dsf = data_source.GetLayer()
    trt_types = []
    for feature in dsf:
        trt_types.append(feature.GetField('TRT'))
    dsf = shapefile.Reader(domains_shapefile)
    dom_shapes = dsf.shapes()
    trt_list = []
    for fault_trace in fault_traces:
        trace_trt_list = []
        for zone_trt, dom_shape in zip(trt_types, dom_shapes):
            dom_poly = Polygon(dom_shape.points)
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
            trt_list.append(default_trt) # Default value 
    return trt_list

def append_xml_header(output_xml,
                      source_model_name):
    """Add information to the nrml which goes above the fault contours

    """
    # Various header info which the user does not control
    output_xml.append("<?xml version='1.0' encoding='utf-8'?>")
    output_xml.append('<nrml xmlns:gml="http://www.opengis.net/gml"' +
                      ' xmlns="http://openquake.org/xmlns/nrml/0.4">')

    # Source model information
    output_xml.append('  <sourceModel name="' + str(source_model_name) + '">')
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
    """Append the fault geometry to the nrml
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
    output_xml.append('          <dip>' + str(dip) + '</dip>')
    output_xml.append('          <upperSeismoDepth>' + str(upper_depth) + '</upperSeismoDepth>')
    output_xml.append('          <lowerSeismoDepth>' + str(lower_depth) + '</lowerSeismoDepth>')
    output_xml.append('      </simpleFaultGeometry>')
    return

def append_incremental_mfd(output_xml, magnitude_scale_rel,
                           rupture_aspect_ratio, rake,
                           min_mag, bin_width, rates,
                           fault_type='simple'):
    """Creates output xml in increment MFD format from 
    given rate parameters
    """
    output_xml.append('      <magScaleRel>' +
                      str(magnitude_scale_rel) + '</magScaleRel>')
    output_xml.append(
        '      <ruptAspectRatio>' + str(rupture_aspect_ratio) + '</ruptAspectRatio>')
    output_xml.append('      <incrementalMFD minMag="' +
                      str(min_mag) + '" binWidth="' + str(bin_width) + '">')
    output_xml.append('          <occurRates>' + ' '.join(str(rt) for rt in rates ) +
                      '</occurRates>')
    output_xml.append('      </incrementalMFD>')
    output_xml.append('      <rake>' + str(rake) + '</rake>')
    if fault_type=='simple':
        output_xml.append('    </simpleFaultSource>')
    elif fault_type=='complex':
        output_xml.append('    </complexFaultSource>')
 
def sliprate2GR_incremental(sliprate, fault_area, b_value,
                            max_mag, min_mag = 0.0, bin_width = 0.1):
    """Converts a sliprate and b-value into a Gutenberg-
    Richter distribution, then converts this to an
    OpenQuake incremental MFD (to facilitate collapsing of rates
    with other MFDs)
    """
    a_value, moment_rate = fault_slip_rate_GR_conversion.slip2GR(sliprate, fault_area,
                                                                 float(b_value), 
                                                                 float(max_mag),
                                                                 M_min=min_mag)

    mfd = TruncatedGRMFD(min_mag, max_mag, bin_width, a_value, b_value)
    mags,rates=zip(*mfd.get_annual_occurrence_rates())
    mags = np.array(mags)
    rates = np.array(rates)
    return mags, rates, moment_rate

def momentrate2YC_incremental(characteristic_mag, b_value,
                              min_mag, max_mag,
                              moment_rate, bin_width):
    """Converts a moment rate and b-value into a Youngs &
    Coppersmith 1985 characteristic distribution, 
    then converts this to an OpenQuake incremental MFD
    (to facilitate collapsing of rates with other MFDs)
    """
    mfd = YoungsCoppersmith1985MFD.from_total_moment_rate(min_mag=0.01,
                                                          b_val=float(b_value),
                                                          char_mag=characteristic_mag, 
                                                          total_moment_rate=moment_rate,
                                                          bin_width=float(bin_width))
    
    mags,rates=zip(*mfd.get_annual_occurrence_rates())
    # calcualate total moment rate and rescale rates if
    # necessary to meet total input rate
    total_moment_rate = 0
    for i in range(len(mags)):
        moment = np.power(10, (1.5*mags[i]+16.05))
        moment = moment/1e7 #Nm
        inc_moment_rate = moment*rates[i]
        total_moment_rate += inc_moment_rate
    moment_error = (total_moment_rate - moment_rate)/moment_rate
    print 'Relative moment rate error', moment_error
    # Rescale rates
    rates = rates/(1+moment_error)
    # Check rates sum as expected
    total_moment_rate = 0
    for i in range(len(mags)):
        moment = np.power(10, (1.5*mags[i]+16.05))
        moment = moment/1e7 #Nm
        inc_moment_rate = moment*rates[i]
        total_moment_rate += inc_moment_rate
    moment_error = (total_moment_rate - moment_rate)/moment_rate
    print 'Final moment rate error',  moment_error

    # Now trim the distribution to just above min_mag
    mags = np.array(mags)
    rates = rates[np.where(mags >= float(min_mag))]
    mags = mags[np.where(mags >= float(min_mag))]
    return mags, rates

def momentrate2MM_incremental(max_mag, moment_rate, bin_width):
   """Converts a moment rate and maximum magmitnude into a
    maximum magnitude distribution as an OpenQuake incremental MFD
    (to facilitate collapsing of rates with other MFDs)
    """
   max_mag = np.around(max_mag, 2)
   mags = np.arange(max_mag - 0.5, max_mag, bin_width)
   print 'mag_values', mags
   moment_values = np.power(10, (1.5*mags+16.05))/1e7
   rate = moment_rate/np.sum(moment_values)
   print 'rate', rate
   
   # check rates sum as expected
   total_moment_rate = 0
   for i in range(len(mags)):
       moment = np.power(10, (1.5*mags[i]+16.05))
       moment = moment/1e7 #Nm
       inc_moment_rate = moment*rate
       total_moment_rate += inc_moment_rate
   moment_error = (total_moment_rate - moment_rate)/moment_rate
   print 'Final moment rate error',  moment_error
   print moment_rate, total_moment_rate
   rates = np.ones(len(mags))*rate
   print 'Total_rate', sum(rates)
   return mags, rates
