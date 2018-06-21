"""Code for building the fault source model with rates from 
collapsed logic tree branches

Jonathan Griffin
Geoscience Australia
April 2017
"""

#module use ~access/modules
#module load pythonlib/basemap

import os, sys
import numpy as np
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_simplefault, \
    shapefile_2_simplefault_CE, shapefile_2_simplefault_MM, shp2nrml
from NSHA2018.source_models.logic_trees import logic_tree
from NSHA2018.source_models.utils.utils import largest_remainder
from NSHA2018.source_models.utils.area_sources import nrml2sourcelist, \
    area2pt_source, weighted_pt_source
from NSHA2018.source_models.utils.pt2fault_distance import read_simplefault_source, \
    pt2fault_distance, write_combined_faults_points, combine_pt_sources, read_pt_source
#from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from openquake.hazardlib.sourcewriter import write_source_model
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from NSHA2018.source_models.utils.scaling import Leonard2014_SCR_extension
from subprocess import call

from openquake.hazardlib.sourceconverter import SourceConverter, \
    area_to_point_sources, SourceGroup

# Basic parameters
shapefile = 'FSM/FSD_simple_faults.shp'
shapefile_faultname_attribute = 'Name'
shapefile_dip_attribute = 'Dip'
shapefile_sliprate_attribute = 'SL_RT_LT'
shapefile_uplift_attribute = 'UP_RT_LT'
source_model_name = 'National_Fault_Source_Model_2018_Collapsed_AUS6_2018'
simple_fault_tectonic_region = None # Define based on neotectonic domains
magnitude_scaling_relation = 'Leonard2014_SCR'
rupture_aspect_ratio = 1.5
upper_depth = 0.001
default_lower_depth = 20.0
a_value = None
#b_region_shapefile =  '../zones/shapefiles/Leonard2008/LEONARD08_NSHA18_MFD.shp'
b_region_shapefile =  '../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
default_b = 1.0#None # Get from Leonard 2008 regions
default_min_mag = 5.5 #4.8
minimum_allowed_min_mag = 4.4 # Protect against problems with rupture mesh spacing
#max_mag = 7.5 #None # Get from scaling
rake = 90
output_dir = source_model_name
#combined_output_dir = 'National_Seismotectonic_Source_Model_2018_Collapsed_NSHA13_2018'
bin_width = 0.1 # Width of MFD bins in magnitude units
domains_shapefile = '../zones/shapefiles/NSHA13_Background/NSHA13_Background_NSHA18.shp'

area_source_model = '../zones/2018_mw/AUS6/input/collapsed/AUS6_collapsed.xml'
#area_source_model = '../zones/2012_mw_ge_4.0/AUS6/input/collapsed/AUS6_collapsed.xml'
#area_source_model = '../zones/2012_mw_ge_4.0/DIMAUS/input/collapsed/DIMAUS_collapsed.xml'
area_source_model_name = area_source_model.split('/')[0].rstrip('.xml')
investigation_time = 50
fault_mesh_spacing = 2 #2 Fault source mesh
rupture_mesh_spacing = 2 #10 # Area source mesh
area_source_discretisation = 15 #20

scalerel = Leonard2014_SCR_extension()
# Define a list of area sources that we won't convert to points, 
# to avoid converting plate boundary region area sources to point
# sources and creating enormous files
area_source_ids = ['NTS', 'TB', 'WBS', 'TAFS', 'WPG', 'PFTB', 'MTB',
                   'BTFZS', 'ADB', 'HP', 'SBS', 'NBOT', 'NBT', 'WM',
                   'WB', 'PFTB-FPB', 'RNBD', 'NBTD', 'BTFZD', 'OSFZD',
                   'MTBD', 'BSSL', 'WRFS', 'CDW', 'WRP', 'NECS', 'NWO',
                   'LGR', 'BS_40_100', 'SRM_40_100', 'BS_100_200', 
                   'BS_200_300', 'SRM_100_200', 'SRM_200_300', 'TMR_40_100',
                   'NT_40_100', 'TMR_100_200', 'NT_100_200', 'BS_300_400',
                   'TMR_200_300', 'NT_200_300', 'TMR_300_400', 'NT_300_400',
                   'BS_400_500', 'BS_500_600', 'BS', 'SWOB', 'WOB', 'SEOB',
                   'NWMB2', 'NWMB1']

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
                                default_trt = 'Non_cratonic')
trt_list = list(set(trts)) # unique trt values
print trt_list
b_values = shp2nrml.b_value_from_region(fault_traces, 
                                        b_region_shapefile, 
                                        default_b = 1.0)
# Output to be appened line by line to this list
output_xml_add = []
output_xml_mb = []
output_xml_geom = []
output_xml_all_methods = []
output_xml_all_methods_inc_cluster = []
output_xmls = [output_xml_add, output_xml_mb, output_xml_geom, output_xml_all_methods, output_xml_all_methods_inc_cluster]
# Append nrml headers
shp2nrml.append_xml_header(output_xml_add, ('%s_additive' % source_model_name))
shp2nrml.append_xml_header(output_xml_mb, ('%s_moment_balanced' % source_model_name))
shp2nrml.append_xml_header(output_xml_geom, ('%s_geom_filtered' % source_model_name))
shp2nrml.append_xml_header(output_xml_all_methods, ('%s_all_methods_collapsed' % source_model_name))
shp2nrml.append_xml_header(output_xml_all_methods_inc_cluster, ('%s_all_methods_collapsed_inc_cluster' % source_model_name))

# Dict of faults where we can apply clustering model based on palaeo-earthquake data
# Rate is present-day ~Mmax rate based on timing of most recent palaeo-earthquake
# Most conservative (i.e. most recent possible) data of palaeoearthqake is used
clustered_fault_rate_dict = {'Cadell Fault':1.399683e-05,
                             'Lake George Scarp': 1.413462e-06,
                             'Hyden Scarp': 1.167942e-05, 
                             'Lake Edgar Fault': 1.40078e-05}



for i, fault_trace in enumerate(fault_traces):
     # Get basic parameters
     dip = dips[i]
     # calculate width based on aspect ratios
     # lengths and widths are in km
     #width = fault_lengths[i]/1.5
     # Calculate width based on Leonard 2014 scaling relations
     length_m = fault_lengths[i]*1000. # convert km to m
     width_m = scalerel.get_width_from_length(length_m, rake)
     width = width_m/1000. # convert back to km
     # Limit width to seismogenic zone
     max_down_dip_width = (default_lower_depth-upper_depth)/np.sin(np.radians(dip))
     print width, max_down_dip_width
     if width > max_down_dip_width:
          width = max_down_dip_width
          lower_depth = default_lower_depth
     elif width < max_down_dip_width:
          lower_depth = width*np.sin(np.radians(dip)) + upper_depth
     print fault_lengths[i], width, dip
     print 'lower_depth', lower_depth
     fault_area = fault_lengths[i]*width #km^2
     #fault_area = fault_lengths[i]*(float(lower_depth)-float(upper_depth))
     sliprate = sliprates[i]
     trt = trts[i]
     faultname = faultnames[i]
     b_value = b_values[i]
     print 'Calculating rates for %s in domain %s' % (faultname, trt)
     # Calculate M_max from scaling relations
     #scalrel = Leonard2014_SCR()
     max_mag = scalerel.get_median_mag(fault_area, float(rake))
     # Round to nearest 0.05 mag unit
     max_mag = np.round((max_mag-0.05), 1) + 0.05
     # Ensure mmax exceeds default Mmin, otherwise we adjust this
     if max_mag < default_min_mag + bin_width:
          min_mag = max_mag - bin_width
     else:
          min_mag = default_min_mag
     print 'Maximum magnitude is %.3f' % max_mag
     if min_mag < minimum_allowed_min_mag:
          print 'Skipping fault %s as not big enough for minimum allowed min_mag of 4.4' % faultname
          continue

     # Append geometry information
     # Use trt type Non_cratonic for Extended
     # only for ground motions, not clustering
     gmm_trt = trt
     if gmm_trt == 'Extended':
          gmm_trt == 'Non_cratonic'
     for output_xml in output_xmls:
          shp2nrml.append_rupture_geometry(output_xml, fault_trace,
                                           dip, i, faultname,
                                           upper_depth, lower_depth,
                                           gmm_trt)
                                           

     # Get truncated Gutenberg-Richter rates
     gr_mags, gr_rates, moment_rate = \
         shp2nrml.sliprate2GR_incremental(sliprate, fault_area,
                                          b_value, max_mag, 
                                          min_mag, bin_width)
         
     gr_mags = np.around(gr_mags, 2)# Rounding to ensure equality of magnitude bins
     gr_add_value, gr_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_Add')
     gr_mb_value, gr_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_MB')
     gr_geom_value, gr_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_Geom')
     gr_add_weight = gr_add_weight[0]
     gr_mb_weight =  gr_mb_weight[0]
     gr_geom_weight = gr_geom_weight[0]

     # Get Youngs and Coppersmith 1985 Characteristic rates
     char_mag = gr_mags[-1] - 0.25 + 0.05 # Adding 0.05 avoids round issues
     print char_mag, b_value, min_mag, max_mag, moment_rate, bin_width
     ce_mags, ce_rates = shp2nrml.momentrate2YC_incremental(char_mag, b_value,
                                                            min_mag, max_mag,
                                                            moment_rate, bin_width)
     ce_mags = np.around(ce_mags, 2)# Rounding to ensure equality of magnitude bins
     ce_add_value, ce_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_Add')
     ce_mb_value, ce_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_MB')
     ce_geom_value, ce_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_Geom')
     ce_add_weight =  ce_add_weight[0]
     ce_mb_weight = ce_mb_weight[0]
     ce_geom_weight = ce_geom_weight[0]

     # Get Maximum Magnitude distirbution
     mm_max_mag = gr_mags[-1] + 0.1 # Avoid rounding issues
     mm_mags, mm_rates = shp2nrml.momentrate2MM_incremental(mm_max_mag, 
                                                            moment_rate,
                                                            bin_width)
     mm_mags = np.around(mm_mags, 2) # Rounding to ensure equality of magnitude bins
     mm_add_value, mm_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_Add')
     mm_mb_value, mm_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_MB')
     mm_geom_value, mm_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_Geom')
     mm_add_weight = mm_add_weight[0]
     mm_mb_weight = mm_mb_weight[0]
     mm_geom_weight = mm_geom_weight[0]
 
     # Get clustered weights for faults with palaeo-earthquake data
     if faultname in clustered_fault_rate_dict:
          cl_value, cl_weight = lt.get_weights('FSM_clustering', trt, branch_value = 'Clustered')
          pois_value, pois_weight = lt.get_weights('FSM_clustering', trt, branch_value = 'Poisson')
          td_value, td_weight = lt.get_weights('FSM_cluster_type', 'Cluster_type', branch_value = 'Time-dependent')
          bimod_value, bimod_weight = lt.get_weights('FSM_cluster_type', 'Cluster_type', branch_value = 'Bimodal')
          # Because we are just calculating mean hazard, bimodal branch collapses to the Poisson branch
          # so weights are combined
          td_weight = td_weight[0]
          pois_weight = pois_weight[0]
          cl_weight = cl_weight[0]
          bimod_weight = bimod_weight[0]
          print td_weight
          print cl_weight
          td_weight = td_weight*cl_weight
          pois_weight = pois_weight + cl_weight*bimod_weight
          print 'td_weight', td_weight
          print 'pois_weight', pois_weight
          cl_mags = mm_mags
          cl_rates = np.ones(len(cl_mags))*(clustered_fault_rate_dict[faultname]/len(cl_mags))
          print cl_rates
          #cl_rates = cl_rates * td_weight
    
     # Calculate collapsed weights
     additive_rates = []
     mb_rates = []
     geom_rates = []
     all_method_rates = []
     all_method_inc_cluster_rates = []
     for mag_bin in gr_mags:
          additive_rate = np.sum(gr_add_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_add_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_add_weight*mm_rates[np.where(mm_mags == mag_bin)])
          additive_rates.append(additive_rate)
          mb_rate = np.sum(gr_mb_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_mb_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_mb_weight*mm_rates[np.where(mm_mags == mag_bin)])
          mb_rates.append(mb_rate)
          geom_rate = np.sum(gr_geom_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_geom_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_geom_weight*mm_rates[np.where(mm_mags == mag_bin)])
          geom_rates.append(geom_rate)
          all_method_rate = np.sum(gr_add_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_add_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_add_weight*mm_rates[np.where(mm_mags == mag_bin)]) + \
              np.sum(gr_mb_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_mb_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_mb_weight*mm_rates[np.where(mm_mags == mag_bin)]) + \
              np.sum(gr_geom_weight*gr_rates[np.where(gr_mags == mag_bin)]) + \
              np.sum(ce_geom_weight*ce_rates[np.where(ce_mags == mag_bin)]) + \
              np.sum(mm_geom_weight*mm_rates[np.where(mm_mags == mag_bin)])
          all_method_rates.append(all_method_rate)
          # Apply clustered rates if we can
          if faultname in clustered_fault_rate_dict:
               all_method_inc_cluster_rate = pois_weight*all_method_rate + \
                   np.sum(td_weight*cl_rates[np.where(cl_mags == mag_bin)])
               all_method_inc_cluster_rates.append(all_method_inc_cluster_rate)
          else:
               all_method_inc_cluster_rates.append(all_method_rate)
               
     # Append collapsed weights to xml
     shp2nrml.append_incremental_mfd(output_xml_add, magnitude_scaling_relation,
                                     rupture_aspect_ratio, rake,
                                     min(gr_mags), bin_width, additive_rates)
     shp2nrml.append_incremental_mfd(output_xml_mb, magnitude_scaling_relation,
                                     rupture_aspect_ratio, rake,
                                     min(gr_mags), bin_width, mb_rates)
     shp2nrml.append_incremental_mfd(output_xml_geom, magnitude_scaling_relation,
                                     rupture_aspect_ratio, rake,
                                     min(gr_mags), bin_width, geom_rates)
     shp2nrml.append_incremental_mfd(output_xml_all_methods, magnitude_scaling_relation,
                                     rupture_aspect_ratio, rake,
                                     min(gr_mags), bin_width, all_method_rates)
     shp2nrml.append_incremental_mfd(output_xml_all_methods_inc_cluster, magnitude_scaling_relation,
                                     rupture_aspect_ratio, rake,
                                     min(gr_mags), bin_width, all_method_inc_cluster_rates)
# Close xml
for output_xml in output_xmls:
     output_xml.append('  </sourceModel>')
     output_xml.append('</nrml>')
# Add newlines
output_xml_add = [oxml + '\n' for oxml in output_xml_add]
output_xml_mb = [oxml + '\n' for oxml in output_xml_mb]
output_xml_geom = [oxml + '\n' for oxml in output_xml_geom]
output_xml_all_methods = [oxml + '\n' for oxml in output_xml_all_methods]
output_xml_all_methods_inc_cluster = [oxml + '\n' for oxml in output_xml_all_methods_inc_cluster]
# Write to file fault models on their own
try:
     os.mkdir(source_model_name)
except:
     pass

f = open(os.path.join(source_model_name, source_model_name + '_additive.xml'),
         'w')
f.writelines(output_xml_add)
f.close()
f = open(os.path.join(source_model_name, source_model_name + '_moment_balance.xml'),
         'w')
f.writelines(output_xml_mb)
f.close()
f = open(os.path.join(source_model_name, source_model_name + '_geom_filtered.xml'),
         'w')
f.writelines(output_xml_geom)
f.close() 
f = open(os.path.join(source_model_name, source_model_name + '_all_methods_collapsed.xml'),
         'w')
f.writelines(output_xml_all_methods)
f.close()
f = open(os.path.join(source_model_name, source_model_name + '_all_methods_collapsed_inc_cluster.xml'),
         'w')
f.writelines(output_xml_all_methods_inc_cluster)
f.close()
#Free memory
del output_xml_add
del output_xml_mb
del output_xml_geom
del output_xml_all_methods
del output_xml_all_methods_inc_cluster

# Now read in the area source model
print 'Reading area source model %s' % area_source_model
area_sources = nrml2sourcelist(area_source_model, 
                          investigation_time=investigation_time, 
                          rupture_mesh_spacing=rupture_mesh_spacing, 
                          width_of_mfd_bin=bin_width,
                          area_source_discretisation=area_source_discretisation)
# Convert area sources to point sources for filtering
print 'Converting to point sources'
area_pt_filename = area_source_model[:-4] + '_pts.xml'
#name = area_source_model.split('/')[-1][:-4] + '_pts'
point_sources, banda_fault_sources, subduction_area_sources = area2pt_source(area_source_model, sources=area_sources,
                                                    filename=area_pt_filename,
                                                    name=source_model_name, return_faults=True,
                                                    exclude_ids = area_source_ids)
pt_source_list = []
for source_group in point_sources:
     for source in source_group:
          pt_source_list.append(source)

# Now  apply weightings diretly to the point source files
# and write one .xml file for each of the methods
# (add, mb, and geom). Note these are summed as we 
# have already collapsed rates for the fault sources
total_add_weight = {}
total_mb_weight = {}
total_geom_weight = {}
for trt in trt_list:
     gr_add_value, gr_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_Add')
     gr_mb_value, gr_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_MB')
     gr_geom_value, gr_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'GR_Geom')
     ce_add_value, ce_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_Add')
     ce_mb_value, ce_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_MB')
     ce_geom_value, ce_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'YC_Geom')
     mm_add_value, mm_add_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_Add')
     mm_mb_value, mm_mb_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_MB')
     mm_geom_value, mm_geom_weight = lt.get_weights('FSM_MFD', trt, branch_value = 'MM_Geom')
     # Get total method weights as these will be applied to zone/pt sources
     total_add_weight[trt] = gr_add_weight[0] + ce_add_weight[0] + mm_add_weight[0]
     total_mb_weight[trt] = gr_mb_weight[0] + ce_mb_weight[0] + mm_mb_weight[0]
     total_geom_weight[trt] = gr_geom_weight[0] + ce_geom_weight[0] + mm_geom_weight[0]
#####################
print 'Not including moment balanced approach for now!'
for trt in trt_list:
     print 'Renormalising weights for other methods'
     partial_weight_sum = total_add_weight[trt]+total_geom_weight[trt]
     total_add_weight[trt] = total_add_weight[trt]*(1/partial_weight_sum)
     total_geom_weight[trt] = total_geom_weight[trt]*(1/partial_weight_sum)
     total_add_weight[trt], total_geom_weight[trt] = largest_remainder([total_add_weight[trt],
                                                                        total_geom_weight[trt]],
                                                                       expected_sum=1,precision=3)
print 'Making weights for subduction regions the same as Non_cratonic'
total_add_weight['Subduction'] = total_add_weight['Non_cratonic']
total_geom_weight['Subduction'] = total_geom_weight['Non_cratonic']
total_mb_weight['Subduction'] = total_mb_weight['Non_cratonic']
##################
additive_pt_sources_filename =  area_source_model[:-4] + '_pts_add_weighted.xml'
model_name = area_source_model.split('/')[-1].rstrip('.xml') + '_additive'
print 'Writing %s' % model_name
additive_pt_sources = weighted_pt_source(pt_source_list, total_add_weight,
                                         model_name, additive_pt_sources_filename, 
                                         nrml_version='04')
mb_pt_sources_filename =  area_source_model[:-4] + '_pts_mb_weighted.xml'
model_name = area_source_model.split('/')[-1].rstrip('.xml') + '_mb'
print 'Writing %s' % model_name
mb_pt_sources = weighted_pt_source(pt_source_list, total_mb_weight,
                                   model_name, mb_pt_sources_filename, 
                                   nrml_version='04')
geom_pt_sources_filename =  area_source_model[:-4] + '_pts_geom_weighted.xml'
model_name = area_source_model.split('/')[-1].rstrip('.xml') + '_geom_filter'
print 'Writing %s' % model_name
geom_pt_sources = weighted_pt_source(pt_source_list, total_geom_weight,
                                     model_name, geom_pt_sources_filename, 
                                     nrml_version='04')
#print 'Exiting here'
#sys.exit()
# Apply geometrical filtering
print 'Applying geometrical filtering - this should be pre-calculated using run_geom_filter.sh!'
#pt2fault_distance(geom_pt_sources, fault_sources, min_distance=5.0,
#                  filename=geom_filtered_pt_source_file,
#                  buffer_distance = 100.,
#                  name=source_model_name)
fsm =  os.path.join(source_model_name, source_model_name + '_geom_filtered.xml')
fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)
geom_filtered_pt_source_file = area_source_model[:-4] + '_pts_geom_weighted_merged_parallel_new_ids.xml'
try:
     geom_filtered_pt_sources = read_pt_source(geom_filtered_pt_source_file)
except IOError:
     msg = 'File %s does not exist, need to run run_geom_filter.sh separately to ' \
         'apply geometrical filter first, then redo_ids.py!' % geom_filtered_pt_source_file
     raise IOError(msg)

outfile = os.path.join(source_model_name, source_model_name + '_geom_filtered_zone.xml')
for fs in banda_fault_sources:
     fault_sources.append(fs)
write_combined_faults_points(geom_filtered_pt_sources, fault_sources,
                             outfile, model_name, area_sources=subduction_area_sources, 
                             nrml_version = '04')

# Apply additive approach                                                                                     
print 'Writing full additive model'
fsm = os.path.join(source_model_name, source_model_name + '_additive.xml')
model_name = source_model_name + '_additive'
outfile =  os.path.join(source_model_name, source_model_name + '_additive_zone.xml')
fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)
for fs in banda_fault_sources:
     fault_sources.append(fs)
write_combined_faults_points(additive_pt_sources, fault_sources,
                             outfile, model_name,  area_sources=subduction_area_sources,
                             nrml_version = '04')

# Merge pt source rates                                                                                       
merged_filename = area_source_model[:-4] + '_pts_geom_add_merged_pts.xml'
model_name = area_source_model.split('/')[-1].rstrip('.xml') + '_add_geom_merged'
combined_pt_sources = combine_pt_sources([additive_pt_sources, geom_filtered_pt_sources],
                                         merged_filename, model_name, nrml_version = '04',
                                         id_location_flag='location')

# Combine merged point sources with merged fault source model
print 'Writing collapsed logic tree seismotectonic model'
fsm = os.path.join(source_model_name, source_model_name + '_all_methods_collapsed_inc_cluster.xml')
model_name = source_model_name + '_' + area_source_model_name + '_collapsed_inc_cluster'
outfile = os.path.join(source_model_name, source_model_name + '_' + \
                            area_source_model_name +'_all_methods_collapsed_inc_cluster.xml')
fault_sources = read_simplefault_source(fsm, rupture_mesh_spacing = fault_mesh_spacing)

# Add Banda fault and area sources extracted earlier
for fs in banda_fault_sources:
     fault_sources.append(fs)
write_combined_faults_points(combined_pt_sources, fault_sources,
                             outfile, model_name,  area_sources=subduction_area_sources,
                             nrml_version = '04')
