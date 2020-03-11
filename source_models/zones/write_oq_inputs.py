from tools.make_nsha_oq_inputs import src_shape2dict, write_oq_sourcefile, \
                                      make_logic_tree, test_beta_curves
from numpy import array, zeros_like, where, unique
from sys import argv
from os import path, mkdir, sep

#print '\n!!! ADD RATE ADJUSTMENT FACTOR !!!\n'

shpfile = argv[1] # input shapefile
outputType = argv[2] # see key below
doSeismotectonic = argv[3] # True = add AU faults; False = ignore AU faults

if doSeismotectonic == 'True':
    doSeismotectonic = True
else:
    doSeismotectonic = False

"""
set outputType:
    0 = best parameters only - use 
    1 = collapsed rates
    2 = multiple file
    3 = multiple beta, best Mmax
    4 = best beta, multiple Mmax

best: best b-value and Mmax (is the same as the highest weighted file in multimod)

collapsed: All variations of b-value and Mmax collapsed into one file

multimod: each combination of b-value and Mmax has its own source file - 15 files in total (should give the same results as above)

mmax_var: collapsed rate file for testing hazard sensitivity by just varying Mmax

bval_var: collapsed rate file for testing hazard sensitivity by varying b-value

"""
    
# get model dictionary from shapefile
model = src_shape2dict(shpfile)
#model = [model[0]]

# split model path
splitpath = shpfile.split(sep)[:-2]

# set beta wts - keep these consistent accross models and zones
beta_wts   = [0.5, 0.2, 0.3] # best, lower (hi b), upper (lo b)

# check that input folder
splitpath.append('input')
modPath = sep.join(splitpath)
  
# check to see if exists
if path.isdir(modPath) == False:
    mkdir(modPath)
        
srcxmls = []

##############################################################################
# build Mmax weights dictionary
##############################################################################

# get Mmax dict from source model weights
ssccsv = '../../shared/seismic_source_model_weights_rounded_p0.4.csv'

lines = open(ssccsv).readlines()

trt = []
mxv = []
mxw = []
for line in lines:
    dat = line.strip().split(',')
    if dat[1] == 'Mmax':
        trt.append(dat[2])
        mxv.append(float(dat[3]))
        mxw.append(float(dat[4]))
        
        # use Proterozoic as generic "Cratonic"
        if dat[2] == 'Proterozoic':
            trt.append('Cratonic')
            mxv.append(float(dat[3]))
            mxw.append(float(dat[4]))
            

trt = array(trt)  
mxv = array(mxv)
mxw = array(mxw)      
unq_trt = unique(trt)

# now make list of dicts
mx_dict = {}
for ut in unq_trt:
    idx = where(trt == ut)[0]
    td = {'trt':ut, 'mx_vals':mxv[idx], 'mx_wts':mxw[idx]}    
    mx_dict[ut] = td
    #print(td

##############################################################################
# use best beta & Mmax
##############################################################################

if outputType == '0':
    
    if doSeismotectonic == True:
        splitpath.append('seismo_best')
        # get output filename
        xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'best_NFSM.xml'
    
    else:
        splitpath.append('best')
        # get output filename
        xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'best.xml'
    
    modPath = sep.join(splitpath)
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
        
    
    # reset beta weights
    tmp_bwts = [1.0, 0.0, 0.0] # only use best beta
    
    #set metadata dict
    meta = {'beta_wts':tmp_bwts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':True, 'mx_idx':-1, 'splitXMLPath': True,
            'doSeisTec':doSeismotectonic} # need to search for best Mmax
    	
    # check to see if exists
    if path.isdir(meta['modelPath']) == False:
        mkdir(meta['modelPath'])
    
    # now write source files
    outxml = write_oq_sourcefile(model, meta, mx_dict)
    srcxmls.append(outxml)
    
    # get branch weight
    branch_wts = [1.0]


##############################################################################
# just write collapsed rates
##############################################################################

elif outputType == '1':
    
    if doSeismotectonic == True:
        splitpath.append('seismo_collapsed')
        # get output filename
        if shpfile.endswith('b_MFD.shp'):
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-6] + 'collapsed_NFSM.xml'
        else:
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed_NFSM.xml'
    
    else:
        splitpath.append('collapsed')
        # get output filename
        if shpfile.endswith('MX.shp'):
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-14] + 'collapsed_mx.xml'
        elif shpfile.endswith('Gridded_b_MFD.shp'):
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-6] + 'collapsed_gridded_b.xml'
        elif shpfile.endswith('b_MFD.shp'):
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-6] + 'collapsed_NFSM.xml'
        else:
            xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.xml'
    
    modPath = sep.join(splitpath)
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
    
    #set metadata dict
    meta = {'beta_wts':beta_wts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':False, 'splitXMLPath': True,
            'doSeisTec':doSeismotectonic}
    	
    # check to see if exists
    if path.isdir(meta['modelPath']) == False:
        mkdir(meta['modelPath'])
    
    # now write source files
    outxml = write_oq_sourcefile(model, meta, mx_dict)
    srcxmls.append(outxml)
    
    # get branch weight
    branch_wts = [1.0]

##############################################################################
# write multiple files
##############################################################################
elif outputType == '2':
    if doSeismotectonic == True:
        splitpath.append('seismo_multimod')
    else:
        splitpath.append('multimod')
    modPath = sep.join(splitpath)
    
    # set Mmax array and weights - cannot equaly weight different TRTs
    mx_wts  = [0.1, 0.2, 0.4, 0.2, 0.1] # using this for default as Mmax probably not major contributor to hazard
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
        
    bSuffix = ['bb', 'bu', 'bl']
    
    branch_wts = []
    
    # loop thru betas    
    for i in range(0, len(beta_wts)):
        
        # set temp beta wt array so full rates are calculated
        tmp_bwts = zeros_like(array(beta_wts))
        tmp_bwts[i] = 1.0
        
        # loop thru mmax - other end will get appropriate TRT
        for j in range(0, len(mx_dict['Archean']['mx_vals'])):
                    
            # get output filename
            xmlfile = '.'.join((path.split(shpfile)[-1].strip('.shp')[:-11],
                                bSuffix[i], 'm'+str(j+1), 
                                'xml'))
            
            #set metadata dict
            meta = {'beta_wts':tmp_bwts, 'modelPath':modPath, 'modelFile':xmlfile, 
                    'multiMods':False, 'one_mx':True, 'mx_idx':j, 'splitXMLPath': True,
                    'doSeisTec':doSeismotectonic}
            	
            # check to see if exists
            if path.isdir(meta['modelPath']) == False:
                mkdir(meta['modelPath'])
            
            # now write source files
            outxml = write_oq_sourcefile(model, meta, mx_dict)
            srcxmls.append(outxml)
            
            # get branch weight
            branch_wts.append(beta_wts[i] * mx_wts[j])

##############################################################################
# vary beta only, mmax = constant
##############################################################################

elif outputType == '3':
    
    splitpath.append('bval_var')
    modPath = sep.join(splitpath)
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
        
    # get output filename
    xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.bvar.xml'
    
    #set metadata dict
    meta = {'beta_wts':beta_wts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':True, 'mx_idx':-1, 'splitXMLPath': True,
            'doSeisTec':doSeismotectonic} # need to search for best Mmax
    	
    # check to see if exists
    if path.isdir(meta['modelPath']) == False:
        mkdir(meta['modelPath'])
    
    # now write source files
    outxml = write_oq_sourcefile(model, meta, mx_dict)
    srcxmls.append(outxml)
    
    # get branch weight
    branch_wts = [1.0]
        
##############################################################################
# vary Mmax only, beta = constant
##############################################################################

elif outputType == '4':
    
    splitpath.append('mmax_var')
    modPath = sep.join(splitpath)
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
        
    # reset beta weights
    tmp_bwts = [1.0, 0.0, 0.0] # only use best beta
    
    # get output filename
    xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.mvar.xml'
    
    #set metadata dict
    meta = {'beta_wts':tmp_bwts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':False, 'splitXMLPath': True, 
            'doSeisTec':doSeismotectonic}
    	
    # check to see if exists
    if path.isdir(meta['modelPath']) == False:
        mkdir(meta['modelPath'])
    
    # now write source files
    outxml = write_oq_sourcefile(model, meta, mx_dict)
    srcxmls.append(outxml)
    
    # get branch weight
    branch_wts = [1.0]
        

##############################################################################
# make logic tree file
##############################################################################
print(branch_wts)
make_logic_tree(srcxmls, branch_wts, meta)

##############################################################################
# plot curves for testing
##############################################################################

#test_beta_curves(model[0], 0.1, meta, mx_dict)

