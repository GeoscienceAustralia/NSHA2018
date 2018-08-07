"""Code for building a source model for plate boundary regions
to the north of Australia
Jonnathan Griffin
Geoscience Australia, November 2017
"""

import os, sys
from glob import glob
import ogr
import numpy as np
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_complexfault, \
    shp2nrml
try:
    import netCDF4
except ImportError:
    print 'NetCDF4 not installed/loaded, be sure to run module load netcdf on NCI'

contour_dir = './shp'
contour_shapefiles = glob(os.path.join(contour_dir, '*.shp'))
rate_dir = './source_rates'
source_model_name = 'Banda_Fault_Sources_NSHA_2018_SimpleFault'
tectonic_region = 'Subduction'
magnitude_scale_rel = 'StrasserInterface'
rupture_aspect_ratio = 1.5
rake = 90 # Fix for normal sources!
default_dip = 30
#min_mag = 6.0 # currently 7.2 based on PTHA curves
bin_width = 0.1
# Store output xml strings here
output_xml = []

# Add xml headers
shp2nrml.append_xml_header(output_xml, source_model_name)

java_limit=104.8 # longitude; western limit of java segment, to clip contours by
source_mmin_dict = {'arutrough':7.0,
                'banda_detachment':7.3,
                'flores':7.5,
#                'hjort':6.5,
#                'macquarienorth':5.5,
                'moresby_trough':6.5,
                'newguinea':7.5,
#                'puysegur':6.7,
                'seram_thrust':7.5,
#                'solomon':7.5,
                'sunda':8.3,
                'tanimbar':7.0,
                'timor':7.0,
                'trobriand':7.5}

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
source_depth_dict = {'arutrough':30,
                    'banda_detachment':20,
                    'flores':30,
                    #                'hjort':6.5,
                    'moresby_trough':30,
                    'newguinea':50,
                    #                'puysegur':6.7,
                    'seram_thrust':40,
                    #                'solomon':7.5,
                    'sunda':50,
                    'tanimbar':25,
                    'timor':25,
                    'trobriand':30}

# Now add each fault source
id_base = 'banda_'
#i = 0
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

    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(shapefile, 0)
    layer = data_source.GetLayer()
#    fault_traces = []
    j=0
    for feature in layer:
        line = feature.GetGeometryRef().GetPoints()
#        fault_traces.append(line)
        pt_list = [list(pts) for pts in line]
        lons=[]
        lats = []
        for pt in pt_list:
            lons.append(pt[0])
            lats.append(pt[1])
        if boundary is not None:
            line = [ (i,j) for (i,j) in zip(lons,lats) if i >= boundary[0] if i <= boundary[1] if j >= boundary[2] if j <= boundary[3]]
        else:
            line = [(i,j) for (i,j) in zip(lons,lats)]
 #       fault_traces.append(line)
        ID = id_base + sourcename + str(j)
        j+=1
        shp2nrml.append_rupture_geometry(output_xml, line, default_dip,
                                         ID, ID, 0, 
                                         source_depth_dict[sourcename],
                                         'Subduction')
        # Now we want to read and sum the rates from the PTHA data
        if sourcename == 'sunda': # only want rates for java segment
            ratefilename = 'EVENT_RATES/output_mw_rate_curve_sunda_java.csv'
        else:
            ratefilename = 'EVENT_RATES/output_mw_rate_curve_%s.csv' % sourcename
        ratedata = np.genfromtxt(ratefilename, delimiter=',', skip_header=1)
        # Get rates above Mmin
        min_mag = source_mmin_dict[sourcename]
        rates = ratedata[:,1]
        magnitudes = ratedata[:,0]
        rates = rates[np.argwhere(magnitudes >= min_mag)].flatten()
        print min_mag
        # Get rid of zero valued rates
        rates = rates[np.argwhere(rates > 0)].flatten()
        rake = source_rake_dict[sourcename]
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
