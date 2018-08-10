'''
Script to make multiple OpenQuake job.ini files and NCI parameter files for running 
one intensity measure at a time
'''

from os import path, sep
from sys import argv

sourcePath = argv[1] # path to model input files
modelType = int(argv[2]) # 0=Final; 1=Background; 2=Regional; 3=Seismotectonic

# get model from path
sourceModel = sourcePath.split(sep)[1]
sourceModelSimple = sourceModel.split('_')[0]

imts = ['PGA', 'SA(0.05)', 'SA(0.1)', 'SA(0.2)', 'SA(0.3)', 'SA(0.5)', 'SA(0.7)', 'SA(1.0)', 'SA(1.5)', 'SA(2.0)', 'SA(4.0)']

print '\nMAKE OPTION FOR mean_hazard_curves= TRUE OR FALSE\n'

# loop thru imts
jobList = []
for i, imt in enumerate(imts):
    # make simple period str
    imtstrp = imt.replace('(', '')
    imtstrp = imtstrp.replace(')', '')
    imtstrp = imtstrp.replace('.', '')
    
    # first, parse job file
    if modelType == 0:
        jobfile = path.join('templates','job_maps_complete.ini')
    else:
        jobfile = path.join('templates','job_maps_templates.ini')
    
    jobtxt = open(jobfile).read()

    # first replace source model
    jobtxt = jobtxt.replace('DOMAINS', sourceModel.upper())
    
    # replace imt string
    jobtxt = jobtxt.replace('"PGA"', '"'+imt+'"')
    
    # replace simple period string
    jobtxt = jobtxt.replace('PGA', imtstrp)
    
    # set source model logic tree 
    jobtxt = jobtxt.replace('Domains_source', sourceModelSimple+'_source')
    
    # make job file
    jobFile = sourcePath + sep + 'job_maps_' + imtstrp + '.ini'
     
    # write file
    f = open(jobFile, 'wb')
    f.write(jobtxt)
    f.close()
    
    ##########################################################################
    # make params file
    ##########################################################################
    
    # do background
    if modelType == 1:
        paramfile = path.join('templates','params_maps_template_background.txt')
    # do regional
    elif modelType == 2:
        paramfile = path.join('templates','params_maps_template_regional.txt')
    # do seismotectonic
    elif modelType == 3:
        paramfile = path.join('templates','params_maps_template_seismotectonic.txt')
    # do final
    elif modelType == 0:
        paramfile = path.join('templates','params_maps_complete_himem.txt')
    
    # first, parse param file
    paramtxt = open(paramfile).read()

    # first replace source model
    paramtxt = paramtxt.replace('Domains_multi_mc', sourceModel)
    
    # replace simple period string
    paramtxt = paramtxt.replace('PGA', imtstrp)
    
    paramtxt = paramtxt.replace('job_maps_PGA.ini', path.split(jobFile)[-1])
    
    # make job file
    if modelType == 0:
        paramFile = sourcePath + sep + 'params_maps_' + imtstrp + '_himem.txt'
    else:
        paramFile = sourcePath + sep + 'params_maps_' + imtstrp + '.txt'
    
    # write file
    f = open(paramFile, 'wb')
    f.write(paramtxt)
    f.close()

