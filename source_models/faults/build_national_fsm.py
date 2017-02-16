"""Build the national scale fault source model from shapefile inputs

Jonathan Griffin
Geoscience Australia, February 2016
"""

import os, sys
from NSHA2018.source_models.faults.shapefile2nrml import shapefile_2_simplefault, \
    shapefile_2_simplefault_CE, shapefile_2_simplefault_MM

shapefile = 'FSM/FSD_simple_faults.shp'

shapefile_faultname_attribute = 'Name'
shapefile_dip_attribute = 'Dip'
shapefile_sliprate_attribute = 'SL_RT_LT'
shapefile_uplift_attribute = 'UP_RT_LT'
source_model_name = 'National_Fault_Source_Model_2018_GR'
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
source_model_file = source_model_name + '_source_model.xml'
source_model_filename = os.path.join(output_dir, source_model_file)
f = open(source_model_filename, 'w')
f.writelines(output_xml_text)
f.close()
source_models = [source_model_file]
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
    newxml += '                    <uncertaintyWeight>'+str(1)+'</uncertaintyWeight>\n'
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
