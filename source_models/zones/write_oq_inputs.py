from tools.make_nsha_oq_inputs import src_shape2dict, write_oq_sourcefile, make_logic_tree
from numpy import array, zeros_like, where, unique
from sys import argv
from os import path, mkdir, sep

shpfile = argv[1] # input shapefile

outputType = argv[2]

'''
set outputType:
    0 = best parameters only - use 
    1 = collapsed rates
    2 = multiple file
    3 = best beta, multiple Mmax
    4 = multiple beta, best Mmax
'''
    
# get model dictionary from shapefile
model = src_shape2dict(shpfile)

# split model path
splitpath = shpfile.split(sep)[:-2]

# set beta wts
beta_wts   = [0.5, 0.2, 0.3]
#beta_wts   = [1.0, 0.0, 0.0]

# set Mmax array and weights
print '\n!!!! temporary override to compare rate collapse method !!!!\n'
mx_wts  = [0.1, 0.2, 0.4, 0.2, 0.1] # for testing

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
        
        # use Proterozoic as generic "cratonic"
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
    td = {'mx_vals':mxv[idx], 'mx_wts':mxw[idx]}    
    mx_dict[ut] = td

##############################################################################
# just write collapsed rates
##############################################################################

if outputType == '1':
    
    splitpath.append('collapsed')
    modPath = sep.join(splitpath)
    
    # check to see if exists
    if path.isdir(modPath) == False:
        mkdir(modPath)
        
    # get output filename
    xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.xml'
    
    #set metadata dict
    meta = {'beta_wts':beta_wts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':False}
    	
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
    splitpath.append('multimod')
    modPath = sep.join(splitpath)
    
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
            
            '''
            # set temp beta wt array
            tmp_mwts = zeros_like(array(mx_wts))
            tmp_mwts[j] = 1.0
            '''
        
            # get output filename
            xmlfile = '.'.join((path.split(shpfile)[-1].strip('.shp')[:-11],
                                bSuffix[i], 'm'+str(j+1), 
                                'xml'))
            
            # make output filename
            #modelFile = path.join('test', xmlfile)
            
            #set metadata dict
            meta = {'beta_wts':tmp_bwts, 'modelPath':modPath, 'modelFile':xmlfile, 
                    'multiMods':False, 'one_mx':True, 'mx_idx':j}
            	
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
        
    branch_wts = []
    
    '''    
    # find and set model Mmax
    mxidx = argmax(array(mx_wts))
    tmp_mx = [mx_vals[mxidx]]
    '''
    
    # get output filename
    xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.bvar.xml'
    
    #set metadata dict
    meta = {'beta_wts':beta_wts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':True, 'mx_idx':-1} # need to search for best Mmax
    	
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
        
    tmp_bwts = [1.0, 0.0, 0.0] # only use best beta
    
    branch_wts = []
            
    # get output filename
    xmlfile = path.split(shpfile)[-1].strip('shp')[:-11] + 'collapsed.mvar.xml'
    
    #set metadata dict
    meta = {'beta_wts':tmp_bwts, 'modelPath':modPath, 'modelFile':xmlfile, 
            'multiMods':False, 'one_mx':False}
    	
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
make_logic_tree(srcxmls, branch_wts, meta)


















        
