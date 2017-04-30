from tools.make_nsha_oq_inputs import src_shape2dict, write_oq_sourcefile
from sys import argv
from os import path, mkdir

shpfile = argv[1]

# get model dictionary from 
model = src_shape2dict(shpfile)

# use best beta, true or false
best_beta = True

# use best Mmax, true or false
best_mx = True


# set Mmax array and weights
mx_vals = [7.3, 7.4, 7.5, 7.6, 7.7] # for testing
mx_wts  = [0.1, 0.2, 0.4, 0.2, 0.1] # for testing

#set metadata dict
meta = {'best_beta':best_beta, 'best_mx':best_mx, 'mx_vals':mx_vals, 'mx_wts':mx_wts, 'modelPath':'test'}
	
# check to see if exists
if path.isdir(meta['modelPath']) == False:
    mkdir(meta['modelPath'])


# now write source files
#write_oq_sourcefile(model, meta)

#modelpath, logicpath, multimods


# now write OQ file
#oqpath = path.join(rootfolder)
#write_oq_sourcefile(model, oqpath, oqpath, multimods, bestcurve)
