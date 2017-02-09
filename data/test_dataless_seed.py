# -*- coding: utf-8 -*-
"""
Created on Thu Feb 09 12:08:12 2017

@author: u56903
"""
from obspy.io.xseed import Parser
from obspy.core import read
from os import path

# read dataless seed volumes
#parser = Parser(path.join('dataless', 'AU.dataless'))
parser = Parser(path.join('dataless', 'AU.seed'))

# parse example mseed data
st = read(path.join('mseed', '200509212246.mseed'))

for tr in st:
    # get record info    
    seedid=tr.get_id()
    channel = tr.stats.channel
    
    start_time = tr.stats.starttime
    
    # check to see if response info exists
    try:
        tr.stats.coordinates = parser.get_coordinates(seedid,start_time)
        
        paz=parser.get_paz(seedid,start_time)
        
        print seedid, paz
    except:
        print 'No info for: ', seedid