"""Build the national scale fault source model from shapefile inputs

Jonathan Griffin
Geoscience Australia, February 2016
"""

import os, sys
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_simplefault, \
    shapefile_2_simplefault_CE, shapefile_2_simplefault_MM
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from subprocess import call
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
b_value = None # Get from Leonard 2008 regions
min_mag = 4.5
max_mag = 7.5 #None # Get from scaling
rake = 90
output_dir = source_model_name
combined_output_dir = 'National_Seismotectonic_Source_Model_2018'

#cmd = "sed -i '/ml\/0.5/c\ml\/0.4' %s" % outfile
#print cmd
###################
# GR fault model
print 'Building Gutenberg-Richter earthquake fault source model'
output_xml_text=shapefile_2_simplefault.nrml_from_shapefile(shapefile,
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
                                                            shapefile_uplift_attribute,
                                                            True)

# Write to file                                                                                                          
try:
    os.mkdir(output_dir)
except:
    pass
source_model_file_GR = source_model_name + '_source_model_GR.xml'
source_model_filename_GR = os.path.join(output_dir, source_model_file_GR)
f = open(source_model_filename_GR, 'w')
f.writelines(output_xml_text)
f.close()

###########################
#CE fault source model
print 'Building Characteristic earthquake fault source model'
output_xml_text=shapefile_2_simplefault_CE.nrml_from_shapefile(shapefile,
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
                                                            True,
                                                            shapefile_uplift_attribute,
                                                            True)

# Write to file                                                                                                          
try:
    os.mkdir(output_dir)
except:
    pass
source_model_file_CE = source_model_name + '_source_model_CE.xml'
source_model_filename_CE = os.path.join(output_dir, source_model_file_CE)
f = open(source_model_filename_CE, 'w')
f.writelines(output_xml_text)
f.close()

#source_models = [source_model_file_GR, source_model_file_CE]
#weights = [0.5, 0.5]

######################################                                                                     
#Maximum magnitude fault source model                                                                      
print 'Building maximum magnitude fault source model'
output_xml_text=shapefile_2_simplefault_MM.nrml_from_shapefile(shapefile,
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
                                                            True,
                                                            shapefile_uplift_attribute,
                                                            True)

# Write to file                                                  
try:
    os.mkdir(output_dir)
except:
    pass
source_model_file_MM = source_model_name + '_source_model_MM.xml'
source_model_filename_MM = os.path.join(output_dir, source_model_file_MM)
f = open(source_model_filename_MM, 'w')
f.writelines(output_xml_text)
f.close()

source_models = [source_model_file_GR, source_model_file_CE, source_model_file_MM]
weights = [0.33, 0.34, 0.33]
######################################  
                
# Now write the source model logic tree file 
#####################################                                                                                   
print 'Writing logic tree file'
newxml = '<?xml version="1.0" encoding="UTF-8"?>\n'
newxml += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
newxml += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
newxml += '    <logicTree logicTreeID="lt1">\n'
newxml += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n'
newxml += '            <logicTreeBranchSet uncertaintyType="sourceModel"\n' \
    '                                branchSetID="bs1">\n\n'

# make branches
for i, branch in enumerate(source_models):
    newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
    newxml += '                    <uncertaintyModel>'+branch+'</uncertaintyModel>\n'
    newxml += '                    <uncertaintyWeight>'+str(weights[i])+'</uncertaintyWeight>\n'
    newxml += '                </logicTreeBranch>\n\n'

newxml += '            </logicTreeBranchSet>\n'
newxml += '        </logicTreeBranchingLevel>\n'
newxml += '    </logicTree>\n'
newxml += '</nrml>'

# write logic tree to file
outxml = os.path.join(output_dir, source_model_name + '_source_model_logic_tree.xml')
f = open(outxml,'w')
f.write(newxml)
f.close()


########################################################
# Now add to NSHA13 zones
######################################################
# Write to file                                                                                                          
try:
    os.mkdir(combined_output_dir)
except:
    pass

seismo_source_models = ['Seismotectonic_model_CE_faults.xml', 'Seismotectonic_model_GR_faults.xml', 
                        'Seismotectonic_model_MM_faults.xml']
nsha_zones = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/zones/NSHA13/NSHA13_collapsed_rates_FF.xml'
parser = nrmlSourceModelParser(nsha_zones)
nsha_source_model = parser.read_file('NSHA13_zones')
parser = nrmlSourceModelParser(source_model_filename_CE)
source_model_CE = parser.read_file('National_Fault_Source_Model_2018_CE')
source_model_list = [nsha_source_model, source_model_CE]

source_list = []
source_model_index = 0
for source_model in source_model_list:
    source_index = 0
    for source in source_model.sources:
        print 'Adding source model ', source.name
        source_list.append(source)
        source_index += 1
    source_model_index += 1
source_model_list[0].sources = source_list
outfile = ("%s/Seismotectonic_source_model_CE_faults.xml") % combined_output_dir
source_model_list[0].serialise_to_nrml(outfile)
#cmd = "sed -i '/ml\/0.5/c\ml\/0.4' %s" % outfile
#print cmd
#call(cmd, shell=False)

##########################################                                                                
# MM faults                                                                                               
#########################################
nsha_zones = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/zones/NSHA13/NSHA13_collapsed_rates_FF.xml'
parser = nrmlSourceModelParser(nsha_zones)
nsha_source_model = parser.read_file('NSHA13_zones')
parser = nrmlSourceModelParser(source_model_filename_MM)
source_model_CE = parser.read_file('National_Fault_Source_Model_2018_MM')
source_model_list = [nsha_source_model, source_model_MM]

source_list = []
source_model_index = 0
for source_model in source_model_list:
    source_index = 0
    for source in source_model.sources:
        print 'Adding source model ', source.name
        source_list.append(source)
        source_index += 1
    source_model_index += 1
source_model_list[0].sources = source_list
outfile = ("%s/Seismotectonic_source_model_MM_faults.xml") % combined_output_dir
source_model_list[0].serialise_to_nrml(outfile)

##########################################
# GR faults
#########################################
nsha_zones = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/zones/NSHA13/NSHA13_collapsed_rates_FF.xml'
parser = nrmlSourceModelParser(nsha_zones)
nsha_source_model = parser.read_file('NSHA13_zones')
parser = nrmlSourceModelParser(source_model_filename_GR)
source_model_GR = parser.read_file('National_Fault_Source_Model_2018_GR')
source_model_list = [nsha_source_model, source_model_GR]

source_list = []
source_model_index = 0
for source_model in source_model_list:
    source_index = 0
    for source in source_model.sources:
        print 'Adding source model ', source.name
        source_list.append(source)
        source_index += 1
    source_model_index += 1
source_model_list[0].sources = source_list
outfile = ("%s/Seismotectonic_source_model_GR_faults.xml") % combined_output_dir
source_model_list[0].serialise_to_nrml(outfile)
#cmd = "sed -i '/nrml\/0.5/c\nrml\/0.4' outfile"
#print cmd
#os.system(cmd)#, shell=False)

source_models = seismo_source_models
weights = [0.5, 0.5]
######################################                                                                                     
# Now write the source model logic tree file 
#####################################                                                                                   
print 'Writing logic tree file'
newxml = '<?xml version="1.0" encoding="UTF-8"?>\n'
newxml += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
newxml += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
newxml += '    <logicTree logicTreeID="lt1">\n'
newxml += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n'
newxml += '            <logicTreeBranchSet uncertaintyType="sourceModel"\n' \
    '                                branchSetID="bs1">\n\n'

# make branches
for i, branch in enumerate(source_models):
    newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
    newxml += '                    <uncertaintyModel>'+branch+'</uncertaintyModel>\n'
    newxml += '                    <uncertaintyWeight>'+str(weights[i])+'</uncertaintyWeight>\n'
    newxml += '                </logicTreeBranch>\n\n'

newxml += '            </logicTreeBranchSet>\n'
newxml += '        </logicTreeBranchingLevel>\n'
newxml += '    </logicTree>\n'
newxml += '</nrml>'

# write logic tree to file
outxml = os.path.join(combined_output_dir, combined_output_dir + '_source_model_logic_tree.xml')
f = open(outxml,'w')
f.write(newxml)
f.close()
