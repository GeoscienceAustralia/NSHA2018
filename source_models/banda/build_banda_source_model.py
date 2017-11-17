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
contour_shapefiles = glob(os.path.join(contour_dir, 'sunda.shp'))
rate_dir = './source_rates'
source_model_name = 'Banda_Sources_NSHA_2018'
tectonic_region = 'banda'
magnitude_scale_rel = 'Leonard2014_Interplate'
rupture_aspect_ratio = 1.5
rake = 90 # Fix for normal sources!
#min_mag = 6.0 # currently 7.2 based on PTHA curves
bin_width = 0.1
# Store output xml strings here
output_xml = []

# Add xml headers
shp2nrml.append_xml_header(output_xml, source_model_name)

java_limit=104.8 # longitude; western limit of java segment, to clip contours by

# Now add each fault source
id_base = 'banda_'
i = 0
for shapefile in contour_shapefiles:
    sourcename = shapefile.split('/')[-1][:-4]
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
                                                     contours)

    # Now we want to read and sum the rates from the PTHA data
    # to get cumulative rates
#    ratefilename = 'all_uniform_slip_earthquake_events_%s.nc' % sourcename
#    ratefile = os.path.join(rate_dir, ratefilename)
#    ratedata = netCDF4.Dataset(ratefile, 'r')
#    magnitudes = ratedata.variables['Mw'][:]
    ratefilename = 'EVENT_RATES/output_mw_rate_curve_%s.csv' % sourcename
    ratedata = np.genfromtxt(ratefilename, delimiter=',', skip_header=1)
    min_mag = ratedata[0][0]
    rates = ratedata[:,1]
    print min_mag
    print rates
#    rates = ratedata.variables['rate_annual'][:]
#    mw_rates = np.vstack([magnitudes, rates])
    

    # These should be sorted by magnitude, but best not to assume
#    mw_rates = mw_rates[mw_rates[:,1].argsort()]
#    mw_rates = np.fliplr(mw_rates) # Descending order for cum-sum calcs
#    mw_list = []
#    cum_rates = []
#    rate_cumsum = 0
#    mag = mw_rates[1][0]

#    j=0
#    for mw in mw_rates[1][:]:
#        if mw == mag:
#            rate_cumsum += mw_rates[0][j]
#        else:
#            cum_rates.append(rate_cumsum)
#            mw_list.append(mag)
#            mag = mw_rates[1][j] # Update to next magnitude increment
#            rate_cumsum += mw_rates[0][j]
#        j+=1

    # Append final rates
#    mw_list.append(mag)
#    cum_rates.append(rate_cumsum)
    # Reverse rates to start at min mag
#    cum_rates = cum_rates[::-1]
    # remove rates = 0
#    cum_rates = list(filter((0.0).__ne__, cum_rates))
#    min_mag = min(mw_list)    

    shp2nrml.append_incremental_mfd(output_xml, magnitude_scale_rel,
                                    rupture_aspect_ratio, rake,
                                    min_mag, bin_width, rates)

# Close xml    
output_xml.append('  </sourceModel>')
output_xml.append('</nrml>')
# Add newlines
output_xml = [oxml + '\n' for oxml in output_xml]

# Write xml file
f = open((source_model_name + '.xml'),'w')
f.writelines(output_xml)
f.close()
