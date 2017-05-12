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

# parse expert elicitation weights
lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.edit.csv')

# get list of source files from area sources with collapsed rates
xmllist = []
xmlpath = []
rootfolder = path.join('..', 'complete_model' )
for root, dirnames, filenames in walk(rootfolder):
    for filename in filenames:
        if filename.endswith('.xml'): #filename.endswith('_collapsed.xml') and path.split(root)[-1] == 'collapsed':
            xmllist.append(filename)
            xmlpath.append(root)


# set up metadata dictionary
modelPath = getcwd() # path where source logit tree is to be saved
meta = {'modelPath': modelPath, 'modelFile':'nsha18_source_model_logic_tree.xml', 
        'splitXMLPath': True} # assume source files in job dir 

# get set weights
src_type, src_wts = lt.get_weights('Source_model', 'Source_type')

# temporarily set regional weights to seimotectonic weights
print '\!!!!REMEMBER TO DELETE SETTING REGIONAL WEIGHT TO SEISMOTECTONIC WEIGHT!!!!\n'
src_wts[2] = src_wts[-1]

# set branch weights
branch_wts = array([])
branch_files = []

# loop throgh source model types
for st, sw in zip(src_type, src_wts):
    print '\n'+st
    src_type_wts = []
    
    orig_st = st
    
    # assume smoothed faults same weight as smoothed seis
    if st == 'Smoothed_faults':
        st = 'Smoothed_seismicity'
        sw = 0.0 # set to zero for now!
        
    # get weights within source type
    models, mod_wts = lt.get_weights('Source_model', st)
    
    # now loop through models within source type
    for mod, mw in zip(models, mod_wts):
        
        # loop through xml files and find matches
        for xl, xp in zip(xmllist, xmlpath):
            
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
                     
                    # add model paths
                    branch_files.append(path.join(xp, xl))
                    
                    #print xl, mod, st, mw
                    
    # re-normalise source type weights if within type neq 1.0
    src_type_wts = array(src_type_wts) / sum(src_type_wts)
    
    # now append source-type weights to branch wts
    branch_wts = hstack((branch_wts, sw * src_type_wts))
                
# do Jono's rounding thingy here...
print branch_wts
# check weights sum to one!
if not sum(branch_wts) == 1.0:
    print '\nWeights do not sum to 1.0!'
    
    # assume testing, so rescale weights
    if sum(branch_wts) < 0.95:
        print '\nAre you testing?  Rescaling weights!!!'
        branch_wts = branch_wts / sum(branch_wts)
        
    
    
# do largest remainder method to make sure numbers 
updated_weights = largest_remainder(branch_wts, expected_sum=1.0, precision=4)

# write source model logic tree
make_logic_tree(branch_files, updated_weights, meta)


