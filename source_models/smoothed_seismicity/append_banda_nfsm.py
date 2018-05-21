from sys import argv
from os import path

ssFilename = argv[1] # filename to append
appendNFSM = argv[2] # True = append; False = ignore

if appendNFSM == 'True':
    appendNFSM = True
else:
    appendNFSM = False
    
# first parse SS file
lines = open(ssFilename).readlines()
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
        
        indo_png_fault_file = path.join('..', 'zones', '2018_mw', 'Java_Banda_PNG', 'input', 'collapsed', 'Java_Banda_PNG_collapsed_faults.xml')
        blines = open(indo_png_fault_file).readlines()[5:-2]
        for bline in blines:
            newxml += bline

        newxml += line
    
    else:
        newxml += line

if appendNFSM == True:
    newFile = ssFilename[:-4]+'_banda_nfsm.xml'
else:
    newFile = ssFilename[:-4]+'_banda.xml'
    	
f = open(newFile, 'wb')
f.write(newxml)
f.close()

#print '\nSkipping Banda Faults\n'
