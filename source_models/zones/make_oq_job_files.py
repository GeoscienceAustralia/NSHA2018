from os import path, sep
from sys import argv

sourcePath = argv[1] # path to model input files

# get model from path
sourceModel = sourcePath.split(sep)[1]
sourceModelSimple = sourceModel.split('_')[0]

imts = ['PGA', 'SA(0.1)', 'SA(0.2)', 'SA(0.3)', 'SA(0.5)', 'SA(0.7)', 'SA(1.0)', 'SA(2.0)', 'SA(4.0)']

print '\nMAKE OPTION FOR mean_hazard_curves= TRUE OR FALSE\n'

# loop thru imts
jobList = []
for imt in imts:
    # make simple period str
    imtstrp = imt.replace('(', '')
    imtstrp = imtstrp.replace(')', '')
    imtstrp = imtstrp.replace('.', '')
    
    # first, parse job file
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
    
    # first, parse param file
    paramfile = path.join('templates','params_maps_template.txt')
    paramtxt = open(paramfile).read()

    # first replace source model
    paramtxt = paramtxt.replace('Domains_multi_mc', sourceModel)
    
    # replace simple period string
    paramtxt = paramtxt.replace('PGA', imtstrp)
    
    paramtxt = paramtxt.replace('job_maps_PGA.ini', path.split(jobFile)[-1])
    
    # make job file
    paramFile = sourcePath + sep + 'params_maps_' + imtstrp + '.txt'
    
    # write file
    f = open(paramFile, 'wb')
    f.write(paramtxt)
    f.close()

