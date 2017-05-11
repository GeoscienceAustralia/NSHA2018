# -*- coding: utf-8 -*-
"""
Created on Thu May 11 17:41:07 2017

@author: u56903
"""
from os import path, getcwd, walk
#from tools.make_nsha_oq_inputs import make_logic_tree
from logic_tree import LogicTree

# parse expert elicitation weights
lt = LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.edit.csv')

# get list of source files from area sources
xmllist = []
rootfolder = path.join('..', 'zones', '2012_mw_ge_4.0' )
for root, dirnames, filenames in walk(rootfolder):
    #for filename in filter(filenames, '.pdf'):
    for filename in filenames:
        if filename.endswith('_collapsed.xml') and path.split(root)[-1] == 'collapsed':
            xmllist.append(filename)


# set up metadata dictionary
modelPath = getcwd()
meta = {'modelPath': modelPath, 'modelFile':'nsha18_source_model_logic_tree.xml'}

# get set weights


# now get model weights
bak_mods, bak_wts = lt.get_weights('Source_model', 'Background')
reg_mods, reg_wts = lt.get_weights('Source_model', 'Regional')
sms_mods, sms_wts = lt.get_weights('Source_model', 'Smoothed_seismicity')
smo_mods, smo_wts = lt.get_weights('Source_model', 'Seismotectonic')
smf_mods, smf_wts = lt.get_weights('Source_model', 'Smoothed_seismicity') # w/ faults



#make_logic_tree(srcxmls, branch_wts, meta)