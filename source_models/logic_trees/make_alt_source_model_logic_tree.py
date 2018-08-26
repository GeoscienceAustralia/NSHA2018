# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:41:07 2017

@author: u56903
"""
from os import path, getcwd, walk
from tools.make_nsha_oq_inputs import make_logic_tree
from logic_tree import LogicTree
from source_models.utils.utils import largest_remainder
from numpy import array, hstack, unique
#from shutil import copyfile
from sys import argv

altweightcsv = argv[1] # e.g. seismic_source_model_weights_rounded_p0.4.fault_test.csv
lt_suffix = argv[2] # for naming

xmllist = []
weighted_smoothing = False # for weighting two adaptive with different Mmin

###############################################################################
# copy regional source models
###############################################################################

relpath = path.join('..', 'zones', '2018_mw')
faultpath = path.join('..', 'faults')
destinationPath = 'final'

#print '\n!!! Using Original Magnitudes !!!\n'
#relpath = path.join('..', 'zones', '2012_mx_ge_4.0')

# copy NSHA13
sourceXML = path.join(relpath, 'NSHA13', 'input', 'collapsed', 'NSHA13_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
sourceXML = path.join(relpath, 'AUS6', 'input', 'collapsed', 'AUS6_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
sourceXML = path.join(relpath, 'DIMAUS', 'input', 'collapsed', 'DIMAUS_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy background source models
###############################################################################

# copy NSHA13_Background
sourceXML = path.join(relpath, 'NSHA13_Background', 'input', 'collapsed', 'NSHA13_BACKGROUND_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP
sourceXML = path.join(relpath, 'ARUP', 'input', 'collapsed', 'ARUP_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP_Background
sourceXML = path.join(relpath, 'ARUP_Background', 'input', 'collapsed', 'ARUP_Background_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Domains
sourceXML = path.join(relpath, 'Domains_multi_mc', 'input', 'collapsed', 'Domains_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Leonard08
sourceXML = path.join(relpath, 'Leonard2008', 'input', 'collapsed', 'Leonard2008_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy SinMcC2016
sourceXML = path.join(relpath, 'SinMcC2016', 'input', 'collapsed', 'SIN_MCC_collapsed.xml')
targetXML = path.join('..', 'complete_model', destinationPath, path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy seismotectonic source models
###############################################################################

#relpath = path.join('..', 'faults')
#relpath = path.join('..', 'faults', 'mx') # for testing catalogue with original magnitudes

# copy NSHA13
#sourceXML = path.join(relpath, 'NSHA13', 'input', 'seismo_collapsed', 'NSHA13_collapsed_NFSM.xml')
sourceXML = path.join(faultpath, 'National_Fault_Source_Model_2018_Collapsed_NSHA13_2018', \
                      'National_Fault_Source_Model_2018_Collapsed_NSHA13_2018__all_methods_collapsed_inc_cluster.xml')
targetXML = path.join('..', 'complete_model', destinationPath, 'NFSM_NSHA13_' + path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
#sourceXML = path.join(relpath, 'AUS6', 'input', 'seismo_collapsed', 'AUS6_collapsed_NFSM.xml')
sourceXML = path.join(faultpath, 'National_Fault_Source_Model_2018_Collapsed_AUS6_2018', \
                      'National_Fault_Source_Model_2018_Collapsed_AUS6_2018__all_methods_collapsed_inc_cluster.xml')
targetXML = path.join('..', 'complete_model', destinationPath, 'NFSM_AUS6_' + path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
#sourceXML = path.join(relpath, 'DIMAUS', 'input', 'seismo_collapsed', 'DIMAUS_collapsed_NFSM.xml')
sourceXML = path.join(faultpath, 'National_Fault_Source_Model_2018_Collapsed_DIMAUS_2018', \
                      'National_Fault_Source_Model_2018_Collapsed_DIMAUS_2018__all_methods_collapsed_inc_cluster.xml')
targetXML = path.join('..', 'complete_model', destinationPath, 'NFSM_DIMAUS_' + path.split(sourceXML)[-1])
#copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy smoothed seismicity source models
###############################################################################

# copy GA adaptive
if weighted_smoothing == True:
    print 'not relevant'

# Use four smoothed models
else:
    sourceXML = path.join('..', 'smoothed_seismicity', 'Cuthbertson2018', 'cuthbertson2018_source_model_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'cuthbertson2018_source_model_banda.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'Hall2007_2018', 'Hall2007_2018_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'Hall2007_2018_banda.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA adaptive
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness', \
                          'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA fixed kernel 
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness', \
                          'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])


###############################################################################
# copy smoothed seismicity source models with faults
###############################################################################

# copy GA adaptive
if weighted_smoothing == True:
    print 'not relevant'

# Use four smoothed models with faults
else:
    sourceXML = path.join('..', 'smoothed_seismicity', 'Cuthbertson2018', 'cuthbertson2018_source_model_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'cuthbertson2018_source_model_banda_nfsm.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    sourceXML = path.join('..', 'smoothed_seismicity', 'Hall2007_2018', 'Hall2007_2018_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'Hall2007_2018_banda_nfsm.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #source_models/smoothed_seismicity/GA_adaptive_smoothing_collapsed_K3_single_corner_completeness/ 
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness', \
                          'GA_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_NFSM_adaptive_smoothing_collapsed_K3_single_corner_completeness_banda_nfsm.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
    #GA fixed kernel 
    sourceXML = path.join('..', 'smoothed_seismicity', 'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness', \
                          'GA_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda_nfsm.xml')
    targetXML = path.join('..', 'complete_model', destinationPath, 'GA_NFSM_fixed_smoothing_50_3_collapsed_single_corner_completeness_banda_nfsm.xml')
    #copyfile(sourceXML, targetXML)
    xmllist.append(path.split(targetXML)[-1])
    
###############################################################################
# parse weights file
###############################################################################
'''
if weighted_smoothing == True:
    lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.ss_split.csv')
else:
'''
lt = LogicTree(altweightcsv)

# set up metadata dictionary
modelPath = getcwd() # path where source logic tree is to be saved
outxml = '_'.join(('nsha18_source_model_logic_tree',lt_suffix+'.xml'))
meta = {'modelPath': modelPath, 'modelFile': outxml, 
        'splitXMLPath': True} # assume source files in job dir 

# get set weights
src_type, src_wts = lt.get_weights('Source_model', 'Source_type')

# temporarily set smoothed seis weights to smoothed+faults
print '\!!!!REMEMBER TO DELETE SETTING REGIONAL WEIGHT TO SEISMOTECTONIC WEIGHT!!!!\n'
#src_wts[0] += src_wts[1]

###############################################################################
# recalibrate source type weights
###############################################################################

#src_wts[0] = 0. # smoothed
#src_wts[1] = 0. # smoothed faults

# rescale source types
src_wts = array(src_wts)/sum(array(src_wts))

###############################################################################
# get weights
###############################################################################

# set branch weights
branch_wts = array([])
branch_xml = []
mod_dict = []

# loop throgh source model types
for st, sw in zip(src_type, src_wts):
    print '\n'+st
    src_type_wts = []
    
    orig_st = st
    
    # assume intra-smoothed fault models have same weight as smoothed seis
    '''
    if st == 'Smoothed_faults':
        st = 'Smoothed_seismicity'
        #sw = 0.0 # set to zero for now!
    '''    
    # get weights within source type
    models, mod_wts = lt.get_weights('Source_model', st)
    print models
    print mod_wts
    
    # now loop through models within source type
    for mod, mw in zip(models, mod_wts):
        
        # loop through xml files and find matches
        for xl in xmllist:
            
            if xl.upper().startswith(mod.upper()):
            
                # do a couple of checks
                #print mod, xl, mw
                if mod == 'NSHA13' and xl.upper().startswith('NSHA13_BACKGROUND'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                elif mod == 'NSHA13' and orig_st == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'DIMAUS' and orig_st == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'DIMAUS' and orig_st == 'Regional' and xl.endswith('NFSM.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'AUS6' and orig_st == 'Seismotectonic' and xl.endswith('collapsed.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'AUS6' and orig_st == 'Regional' and xl.endswith('NFSM.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'NSHA13' and orig_st == 'Regional' and xl.endswith('NFSM.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_NFSM_adaptive' and orig_st == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_adaptive' and orig_st == 'Smoothed_faults' and xl.endswith('banda.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_NFSM_fixed' and orig_st == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_fixed' and orig_st == 'Smoothed_faults':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'Cuthbertson' and orig_st == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                elif mod == 'Cuthbertson' and orig_st == 'Smoothed_faults' and xl.endswith('banda.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'Hall' and orig_st == 'Smoothed_seismicity' and xl.endswith('nfsm.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                elif mod == 'Hall' and orig_st == 'Smoothed_faults' and xl.endswith('banda.xml'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                # else, add file to list
                else:
                    # multiply ARUP models by 0.5
                    if mod.startswith('ARUP'):
                       mod_wt = 0.5
                       print '    Modifying ARUP model weight'
                    else:
                       mod_wt = 1.0
                    
                    # append weights within source type
                    src_type_wts.append(mod_wt * mw)
                    
                    # append branch file
                    branch_xml.append(xl)
                     
                    #print xl, mod, st, mw
                    
                    # get models actually added
                    mod_dict.append({'xml':xl, 'model':mod, 'model_wt':mw, 'src_type':st, 'src_wt':sw, 'cml_wt':mw*sw})
                    
    # re-normalise source type weights if within type neq 1.0
    src_type_wts = array(src_type_wts) / sum(src_type_wts)
    
    # now append source-type weights to branch wts
    branch_wts = hstack((branch_wts, sw * array(src_type_wts)))
                
#print branch_wts

# check weights sum to one!
if not sum(branch_wts) == 1.0:
    print '\nWeights do not sum to 1.0!:',sum(branch_wts)
    
    # assume testing, so rescale weights
    if sum(branch_wts) < 0.95:
        print branch_wts
        print '\nAre you testing?  Rescaling weights!!!\n'
        #branch_wts = branch_wts / sum(branch_wts)
        
       
# do largest remainder method to make sure numlers 
updated_weights = largest_remainder(branch_wts, expected_sum=1.0, precision=3)
#updated_weights = branch_wts

# write source model logic tree
make_logic_tree(branch_xml, updated_weights, meta)

