# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:41:07 2017

@author: u56903
"""
from os import path, getcwd, walk
from tools.make_nsha_oq_inputs import make_logic_tree
from logic_tree import LogicTree
from source_models.utils.utils import largest_remainder
from numpy import array, hstack
from shutil import copyfile

xmllist = []

###############################################################################
# copy regional source models
###############################################################################

relpath = path.join('..', 'zones', '2012_mw_ge_4.0')

# copy NSHA13
sourceXML = path.join(relpath, 'NSHA13', 'input', 'collapsed', 'NSHA13_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
sourceXML = path.join(relpath, 'AUS6', 'input', 'collapsed', 'AUS6_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
sourceXML = path.join(relpath, 'DIMAUS', 'input', 'collapsed', 'DIMAUS_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy background source models
###############################################################################

# copy NSHA13_Background
sourceXML = path.join(relpath, 'NSHA13_Background', 'input', 'collapsed', 'NSHA13_BACKGROUND_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP
sourceXML = path.join(relpath, 'ARUP', 'input', 'collapsed', 'ARUP_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy ARUP_Background
sourceXML = path.join(relpath, 'ARUP_Background', 'input', 'collapsed', 'ARUP_Background_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Domains
sourceXML = path.join(relpath, 'Domains', 'input', 'collapsed', 'Domains_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy Leonard08
sourceXML = path.join(relpath, 'Leonard08', 'input', 'collapsed', 'Leonard08_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy SinMcC2016
sourceXML = path.join(relpath, 'SinMcC2016', 'input', 'collapsed', 'SinMcC2016_collapsed.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy seismotectonic source models
###############################################################################

relpath = path.join('..', 'faults')

# copy NSHA13
sourceXML = path.join(relpath, 'NFSM_NSHA13_collapsed_additive_w1', 'NFSM_NSHA13_collapsed_additve_w1.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy AUS6
sourceXML = path.join(relpath, 'NFSM_AUS6_collapsed_additive_w1', 'NFSM_AUS6_collapsed_additive_w1.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy DIMAUS
sourceXML = path.join(relpath, 'NFSM_DIMAUS_collapsed_additive_w1', 'NFSM_DIMAUS_collapsed_additive_w1.xml')
targetXML = path.join('..', 'complete_model', path.split(sourceXML)[-1])
copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

###############################################################################
# copy smoothed seismicity source models
###############################################################################
'''
# copy GA fixed
sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_fixed_smoothing_collapsed/source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert_v1.xml'
targetXML = path.join('..', 'complete_model', 'GA_fixed_smoothing_full_uncert.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])
'''

# copy GA adaptive
#sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_adaptive_smoothing_collapsed_K4_mmin3p0/Australia_Adaptive_K4_merged_inc_b_mmax_uncert.xml'
#targetXML = path.join('..', 'complete_model', 'GA_adaptive_smoothing_full_uncert.xml')
sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/Australia_Adaptive_K4_merged_bestb_mmin3.0_w1.xml' # best b
targetXML = path.join('..', 'complete_model', 'GA_M3_adaptive_smoothing_best_b.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'source_model_Australia_Adaptive_K3_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# copy GA adaptive
#sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_adaptive_smoothing_collapsed_K4_mmin3p0/Australia_Adaptive_K4_merged_inc_b_mmax_uncert.xml'
#targetXML = path.join('..', 'complete_model', 'GA_adaptive_smoothing_full_uncert.xml')
sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_adaptive_smoothing_collapsed_K4_mmin4p0_bestb/Australia_Adaptive_K4_merged_bestb_mmin4.0_w1.xml' # best b
targetXML = path.join('..', 'complete_model', 'GA_M3_adaptive_smoothing_best_b.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'source_model_Australia_Adaptive_K3_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])


###############################################################################
# copy smoothed seismicity source models with faults
###############################################################################
'''
# copy GA fixed
sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_fixed_smoothing_collapsed/source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert_v1.xml'
targetXML = path.join('..', 'complete_model', 'GA_NFSM_fixed_smoothing_full_uncert.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'NFSM_source_model_smoothed_frankel_50_3_mmin_3.0_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])
'''

# M3
# copy GA adaptive
#sourceXML = '/short/w84/NSHA18/sandpit/jdg547/NSHA2018/source_models/smoothed_seismicity/GA_adaptive_smoothing_collapsed_K4_mmin3p0_faults/Australia_Adaptive_K4_merged_inc_b_mmax_uncert_faults_additive.xml'
#targetXML = path.join('..', 'complete_model', 'GA_NFSM_adaptive_smoothing_full_uncert.xml')
sourceXML = '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed/GA_adaptive_smoothing_best/Adaptive_best_b_Joined_K4_faults/maps_k4_joined_20170522_225102_Adaptive_best_b_Joined_K4_faults_tia547/Australia_Adaptive_K4_merged_bestb_mmin3.0_w1_faults.xml'
targetXML = path.join('..', 'complete_model', 'GA_NFSM_M3_adaptive_smoothing_best_b.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'NFSM_source_model_Australia_Adaptive_K3_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])

# M4
sourceXML = '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed/GA_adaptive_smoothing_best/Adaptive_best_b_Joined_K4_faults/maps_k4_joined_20170522_225102_Adaptive_best_b_Joined_K4_faults_tia547/Australia_Adaptive_K4_merged_bestb_mmin4.0_w1_faults.xml'
targetXML = path.join('..', 'complete_model', 'GA_NFSM_M4_adaptive_smoothing_best_b.xml')
try:
    copyfile(sourceXML, targetXML)
# FOR TESTING ONLY
except:
    sourceXML = path.join('..', 'testing', 'NFSM_source_model_Australia_Adaptive_K3_merged_inc_b_mmax_uncert_v1.xml')
    copyfile(sourceXML, targetXML)
xmllist.append(path.split(targetXML)[-1])


###############################################################################
# parse weights file
###############################################################################

#lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.edit.csv')
lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.ss_split.csv')

# set up metadata dictionary
modelPath = getcwd() # path where source logit tree is to be saved
meta = {'modelPath': modelPath, 'modelFile':'nsha18_source_model_logic_tree.xml', 
        'splitXMLPath': True} # assume source files in job dir 

# get set weights
src_type, src_wts = lt.get_weights('Source_model', 'Source_type')

# temporarily set smoothed seis weights to smoothed+faults
print '\!!!!REMEMBER TO DELETE SETTING REGIONAL WEIGHT TO SEISMOTECTONIC WEIGHT!!!!\n'
#src_wts[0] += src_wts[1]

###############################################################################
# recalibrate source type weights
###############################################################################

src_wts[0] = 0. # smoothed
src_wts[1] = 0. # smoothed faults

# rescale source types
src_wts = array(src_wts)/sum(array(src_wts))

###############################################################################
# get weights
###############################################################################

# set branch weights
branch_wts = array([])
branch_xml = []

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
    
    # now loop through models within source type
    for mod, mw in zip(models, mod_wts):
        
        # loop through xml files and find matches
        for xl in xmllist:
            
            if xl.upper().startswith(mod.upper()):
            	
                # do a couple of checks
                if mod == 'NSHA13' and xl.upper().startswith('NSHA13_BACKGROUND'):
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                elif mod == 'DIMAUS' and orig_st == 'Seismotectonic':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'NSHA13' and orig_st == 'Seismotectonic':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'AUS6' and orig_st == 'Seismotectonic':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_adaptive' and orig_st == 'Smoothed_faults':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                    
                elif mod == 'GA_fixed' and orig_st == 'Smoothed_faults':
                    print 'Not adding '+xl+' to '+orig_st+' set'
                
                # else, add file to list
                else:
                    # multiply ARUP models by 0.5
                    if mod.startswith('ARUP'):
                       mod_wt = 0.5
                       print 'Modifying ARUP model weight'
                    else:
                       mod_wt = 1.0
                    
                    # append weights within source type
                    src_type_wts.append(mod_wt * mw)
                    
                    # append branch file
                    branch_xml.append(xl)
                     
                    print xl, mod, st, mw
                    
    # re-normalise source type weights if within type neq 1.0
    src_type_wts = array(src_type_wts) / sum(src_type_wts)
    
    # now append source-type weights to branch wts
    branch_wts = hstack((branch_wts, sw * array(src_type_wts)))
                
print branch_wts

# check weights sum to one!
if not sum(branch_wts) == 1.0:
    print '\nWeights do not sum to 1.0!'
    
    # assume testing, so rescale weights
    if sum(branch_wts) < 0.95:
        print '\nAre you testing?  Rescaling weights!!!\n'
        branch_wts = branch_wts / sum(branch_wts)
        
       
# do largest remainder method to make sure numbers 
updated_weights = largest_remainder(branch_wts, expected_sum=1.0, precision=3)
#updated_weights = branch_wts

# write source model logic tree
make_logic_tree(branch_xml, updated_weights, meta)

