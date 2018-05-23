"""Code for building a source model for plate boundary regions
to the north of Australia
Jonnathan Griffin
Geoscience Australia, November 2017
"""

import os, sys
from glob import glob
import numpy as np
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_complexfault, \
    shp2nrml
try:
    import netCDF4
except ImportError:
    print 'NetCDF4 not installed/loaded, be sure to run module load netcdf on NCI'

contour_dir = './shp'
contour_shapefiles = glob(os.path.join(contour_dir, '*.shp'))
#rate_dir = './EVENT_RATES'
source_model_name = 'Banda_Fault_Sources_NSHA_2018'
tectonic_region = 'Subduction'
magnitude_scale_rel = 'StrasserInterface'
rupture_aspect_ratio = 1.5
rake = 90 # Fix for normal sources!
#min_mag = 6.0 # currently 7.2 based on PTHA curves
bin_width = 0.1
# Store output xml strings here
output_xml = []

# Add xml headers
shp2nrml.append_xml_header(output_xml, source_model_name)

java_limit=104.8 # longitude; western limit of java segment, to clip contours by
source_mmin_dict = {'arutrough':6.5,
                'banda_detachment':7.0,
                'flores':7.0,
#                'hjort':6.5,
#                'macquarienorth':5.5,
                'moresby_trough':6.0,
                'newguinea':7.0,
#                'puysegur':6.7,
                'seram_thrust':7.0,
#                'solomon':7.5,
                'sunda':8.0,
                'tanimbar':6.5,
                'timor':6.5,
                'trobriand':6.5}

source_rake_dict = {'arutrough':-90,
                    'banda_detachment':-90,
                    'flores':90,
                    #                'hjort':6.5,
                    'moresby_trough':90,
                    'newguinea':90,
                    #                'puysegur':6.7,
                    'seram_thrust':90,
                    #                'solomon':7.5,
                    'sunda':90,
                    'tanimbar':90,
                    'timor':90,
                    'trobriand':-90}

# Now add each fault source
id_base = 'banda_'
i = 0
for shapefile in contour_shapefiles:
    sourcename = shapefile.split('/')[-1][:-4]
    if sourcename not in source_mmin_dict.keys():
        print 'Skipping %s source zone' % sourcename
        continue
    print 'Adding source %s' % sourcename
    if sourcename == 'sunda': 
        boundary = [java_limit, 360, -90, 90]
    else:
        boundary = None
    contours = shapefile_2_complexfault.parse_line_shapefile(
        shapefile, 'level', boundary=boundary)
    ID = id_base + str(i)
    i+=1
    shapefile_2_complexfault.append_fault_source_header(output_xml,
                                                        ID,
                                                        sourcename,
                                                        tectonic_region)
    shapefile_2_complexfault.append_rupture_geometry(output_xml, 
                                                     contours, dps=4)

    # Now we want to read and sum the rates from the PTHA data
    if sourcename == 'sunda': # only want rates for java segment
        ratefilename = 'EVENT_RATES/output_mw_rate_curve_sunda_java.csv'
    else:
        ratefilename = 'EVENT_RATES/output_mw_rate_curve_%s.csv' % sourcename
    ratedata = np.genfromtxt(ratefilename, delimiter=',', skip_header=1)
    # Get rates above Mmin
    min_mag = source_mmin_dict[sourcename]
    cum_rates = ratedata[:,1]
    magnitudes = ratedata[:,0]
    cum_rates = cum_rates[np.argwhere(magnitudes >= min_mag)].flatten()
    print min_mag
    # Get rid of zero valued rates
    cum_rates = cum_rates[np.argwhere(cum_rates > 0)].flatten()
    # Convert rates to increments, as at present they are cumulative, 
    # i.e. prob Mw >= m
    cum_rates = cum_rates[::-1]
    rates = []
    for j, rate in enumerate(cum_rates):
        if j==0:
            rates.append(rate)
            tmp_rate = rate
        else:
            inc_rate = rate - tmp_rate
            rates.append(inc_rate)
            tmp_rate = rate
    rates = np.array(rates)
    rates = rates[::-1]
    rake = source_rake_dict[sourcename]
    # Adjust min_mag for Openquake MFD distribution
    min_mag += bin_width/2.
    shp2nrml.append_incremental_mfd(output_xml, magnitude_scale_rel,
                                    rupture_aspect_ratio, rake,
                                    min_mag, bin_width, rates,
                                    fault_type='complex')

# Close xml    
output_xml.append('  </sourceModel>')
output_xml.append('</nrml>')
# Add newlines
output_xml = [oxml + '\n' for oxml in output_xml]

# Write xml file
f = open((source_model_name + '.xml'),'w')
f.writelines(output_xml)
f.close()
