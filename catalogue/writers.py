# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:42:50 2017

@author: u56903
"""

def ggcat2ascii(ggcat_dict, outfile):    
    '''
    # takes event dictionary format as parsed by catalogues.parsers.parse_ggcat
    
    # exports to outfile
    '''
    # make file for writing
    cattxt = ''
    
    for ev in ggcat_dict:
        # set string constant width
        datestr = '{0.year:4d} {0.month:02d} {0.day:02d} {0.hour:02d}{0.minute:02d}'.format(ev['datetime'])              
        
        # make line        
        line = ' '.join((datestr, str("%0.3f" % ev['lon']), str("%0.3f" % ev['lat']), \
                         str("%0.1f" % ev['dep']), str("%0.1f" % ev['prefmag']), \
                         ev['prefmagtype'], ev['auth']))
                       
        cattxt = cattxt + line + '\n'
                
    #f = open(path.join(outfolder, outfile), 'wb')
    f = open(outfile, 'wb')
    f.write(cattxt)
    f.close()