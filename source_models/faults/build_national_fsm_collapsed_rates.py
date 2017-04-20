"""Code for building the fault source model with rates from 
collapsed logic tree branches

Jonathan Griffin
Geoscience Australia
April 2017
"""

import os, sys
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_simplefault, \
    shapefile_2_simplefault_CE, shapefile_2_simplefault_MM, shp2nrml
from NSHA2018.source_models.logic_trees import logic_tree
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from subprocess import call

# Basic parameters
shapefile = 'FSM/FSD_simple_faults.shp'
shapefile_faultname_attribute = 'Name'
shapefile_dip_attribute = 'Dip'
shapefile_sliprate_attribute = 'SL_RT_LT'
shapefile_uplift_attribute = 'UP_RT_LT'
source_model_name = 'National_Fault_Source_Model_2018'
simple_fault_tectonic_region = None # Define based on neotectonic domains
magnitude_scaling_relation = 'Leonard2014_SCR'
rupture_aspect_ratio = 1
upper_depth = 0.001
lower_depth = 20.0
a_value = None
b_region_shapefile =  '../zones/shapefiles/Leonard2008/LEONARD08_NSHA18_MFD.shp'
default_b = 1.0#None # Get from Leonard 2008 regions
min_mag = 4.5
#max_mag = 7.5 #None # Get from scaling
rake = 90
output_dir = source_model_name
combined_output_dir = 'National_Seismotectonic_Source_Model_2018'
bin_width = 0.1 # Width of MFD bins in magnitude units
domains_shapefile = '../zones/shapefiles/Domains/DOMAINS_NSHA18.shp'

# Get logic tree information
lt = logic_tree.LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')

# Get basic information from shapefile
fault_traces, faultnames, dips, sliprates, fault_lengths = \
    shp2nrml.parse_line_shapefile(shapefile, 
                                  shapefile_faultname_attribute,
                                  shapefile_dip_attribute, 
                                  shapefile_sliprate_attribute,
                                  shapefile_uplift_attribute=shapefile_uplift_attribute,
                                  slip_units = 'm/ma')

# Get b-value and trt from domains
trts = shp2nrml.trt_from_domains(fault_traces, domains_shapefile,
                                default_trt = 'Non-cratonic')
b_values = shp2nrml.b_value_from_region(fault_traces, 
                                        b_region_shapefile, 
                                        default_b = 1.0)

for i, fault_trace in enumerate(fault_traces):
     fault_area = fault_lengths[i]*(float(lower_depth)-float(upper_depth))
     sliprate = sliprates[i]
     trt = trts[i]
     faultname = faultnames[i]
     b_value = b_values[i]
     print 'Calculating rates for %s in domain %s' % (faultname, trt)
     # Calculate M_max from scaling relations
     scalrel = Leonard2014_SCR()
     max_mag = scalrel.get_median_mag(fault_area, float(rake))
     print 'Maximum magnitude is %.3f' % max_mag

     # Get b-value and trt from domains
#     trt = shp2nrml.trt_from_domains

     # Get truncated Gutenberg-Richter rates
     gr_mags, gr_rates, moment_rate = \
         shp2nrml.sliprate2GR_incremental(sliprate, fault_area,
                                          b_value, max_mag, 
                                          min_mag, bin_width)
     print gr_mags
     print gr_rates
     gr_value, gr_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_Add')
     print gr_weight
     # Get Youngs and Coppersmith 1985 Characteristic rates
     char_mag = gr_mags[-1] - 0.25 + 0.05 # Adding 0.05 avoids round issues
     yc_mags, yc_rates = shp2nrml.momentrate2YC_incremental(char_mag, b_value,
                                                            min_mag, max_mag,
                                                            moment_rate, bin_width)
     print yc_mags
     print yc_rates
     

     # Get Maximum Magnitude distirbution
     mm_max_mag = gr_mags[-1] + 0.1 # Avoid rounding issues
     mm_mags, mm_rates = shp2nrml.momentrate2MM_incremental(mm_max_mag, 
                                                            moment_rate,
                                                            bin_width)
     print mm_mags, mm_rates
     sys.exit()
                                                          


