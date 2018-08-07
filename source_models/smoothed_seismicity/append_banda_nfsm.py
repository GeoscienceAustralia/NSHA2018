'''
Appends finite fault and banda sea sources to smoothed seismicity models

Takes a tar.gz file and unzips to generate output xml
'''

from sys import argv
from os import path, system

xmlFile = argv[1] # xml file to append sources - assume already untarred
appendNFSM = argv[2] # True = append; False = ignore

if appendNFSM == 'True':
    appendNFSM = True
else:
    appendNFSM = False
    
# first de-compress tar.gz file
'''
print 'Extracting tar.gz file ...'
tarFolder = path.split(ssTarFile)[0]
system('cd '+tarFolder)
system('tar -xvf '+ path.split(ssTarFile)[-1])
system('cd ..')
'''

# move xml file to whence it came
#xmlFile = ssTarFile.strip('.tar.gz')
#system(' '.join(('mv', path.split(xmlFile)[-1], path.split(xmlFile)[0]+path.sep)))

# now parse SS file
lines = open(xmlFile).readlines()
newxml = ''

for line in lines:
    # inset additional sources
    if line.strip().startswith('</sourceModel>'):
        
        ######################################################################
        # add Australian fault-source model
        ######################################################################
        if appendNFSM == True:
 
            aust_fault_file = path.join('..', 'faults', 'NFSM', \
                                        'National_Fault_Source_Model_2018_Collapsed_all_methods_collapsed_inc_cluster.xml')

            flines = open(aust_fault_file).readlines()[3:-2]
            for fline in flines:
                newxml += '    ' + fline
        
        ######################################################################
        # add indoneasia-png source model - do this for all
        ######################################################################
        
        #indo_png_fault_file = path.join('..', 'zones', '2018_mw', 'Java_Banda_PNG', 'input', 'collapsed', 'Java_Banda_PNG_collapsed_faults.xml')
        indo_png_fault_file = path.join('..', 'zones', '2018_mw', 'Java_Banda_PNG', 'input', 'collapsed', 'Java_Banda_PNG_collapsed.xml') # file above is replaced (includes faults)
        blines = open(indo_png_fault_file).readlines()[4:-2]
        for bline in blines:
            newxml += bline

        newxml += line
    
    else:
        newxml += line

if appendNFSM == True:
    newFile = xmlFile[:-4]+'_banda_nfsm.xml'
else:
    newFile = xmlFile[:-4]+'_banda.xml'
    	
f = open(newFile, 'wb')
f.write(newxml)
f.close()

#print '\nSkipping Banda Faults\n'
