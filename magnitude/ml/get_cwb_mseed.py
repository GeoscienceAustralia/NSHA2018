
# coding: utf-8

# In[7]:

from catalogue.parsers import parse_ggcat

#!/opt/anaconda-2.3.0/bin/python2
"""#!/opt/antelope/5.4/bin/python"""
#%matplotlib inline
import numpy as np
import datetime
from os import path, makedirs
from obspy.core import utcdatetime, event
from obspy.core.event import Catalog, Event, Magnitude, Origin, StationMagnitude
#from obspy.neic.client import Client
from obspy.clients.neic.client import Client
from obspy.io.xseed import Parser
#from obspy.core.util import gps2DistAzimuth
from obspy.geodetics import gps2dist_azimuth as gps2DistAzimuth
#from obspy.taup import TauPyModel
#from obspy.signal import invsim as inv
#from obspy.io.xseed.utils import SEEDParserException
# plot all traces into one pdf file
#from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib import pyplot

# we will use dataless seed from IRIS to get station information
parser = Parser("../../data/AU.dataless")

# travel-time model will be iasp91 but could be any
from obspy import taup
vel_model = taup.TauPyModel(model="iasp91")

#from obspy.taup.TauPyModel import get_travel_times
# local modules
#import local_magnitude

##########################################################################
# set constants and funcs
##########################################################################

r_earth = 6371.

def sind(x): return np.sin(x / 180. * np.pi)
def cosd(x): return np.cos(x / 180. * np.pi)
def tand(x): return np.tan(x / 180. * np.pi)
def arcsind(x): return np.arcsin(x) / np.pi * 180
def arccosd(x): return np.arccos(x) / np.pi * 180
def arctand(x): return np.arctan(x) / np.pi * 180
def gps2DistDegree(lat1, lon1, lat2, lon2):
    return arccosd(sind(lat1) * sind(lat2) +
                   cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))

##########################################################################
# here we define how we measure peak to peak amplitude and period
##########################################################################
# 
def max_p2t(data, delta):
     """
     Function to find the maximum peak-to-trough amplitude and period of this \
     amplitude.

     :type data: ndarray
     :param data: waveform trace to find the peak-to-trough in.
     :type delta: float
     :param delta: Sampling interval in seconds

     :returns: tuple of (amplitude, period, time) with amplitude in the same \
         scale as given in the input data, and period in seconds, and time in \
         seconds from the start of the data window.
     """
     turning_points = []  # A list of tuples of (amplitude, sample)
     for i in range(1, len(data) - 1):
         if (data[i] < data[i-1] and data[i] < data[i+1]) or\
            (data[i] > data[i-1] and data[i] > data[i+1]):
             turning_points.append((data[i], i))
     if len(turning_points) >= 1:
         amplitudes = np.empty([len(turning_points)-1],)
         half_periods = np.empty([len(turning_points)-1],)
     else:
         print('Turning points has length: '+str(len(turning_points)) +
               ' data have length: '+str(len(data)))
         return (0.0, 0.0, 0.0)
     for i in range(1, len(turning_points)):
         half_periods[i-1] = (delta * (turning_points[i][1] -
                                       turning_points[i-1][1]))
         amplitudes[i-1] = np.abs(turning_points[i][0]-turning_points[i-1][0])
     amplitude = np.max(amplitudes)
     period = 2 * half_periods[np.argmax(amplitudes)]

     return (amplitude, period, delta*turning_points[np.argmax(amplitudes)][1])



# initialize the cwb port
client=Client(host='10.7.161.60',port=2061,debug=False, timeout=60)
eq=[]

# Instantiate catalogue object
catalogue = Catalog()
# Build Quakeml Event object

""" 
JG's original method
event_num = 0
with open("../../eq2.txt",'r') as cat:
    for line in cat:
        
        eq=line.split(',')
        yr=int(eq[4])
        mon=int(eq[5])
        day=int(eq[6])
        hr=int(eq[7])
        mn=int(eq[8])
        sec=float(eq[9])
        code=eq[10]
        corr=float(eq[11])
        lon=float(eq[12])
        lat=-1*float(eq[13]) #Latitudes in Gary Gibson catalogue given in southern hemisphere coordindates
        dep=float(eq[14])
        mag_type = eq[16]
        mag_pref = float(eq[17])
        mag_source = eq[31].split(' ')[-1]
        
        start_time=utcdatetime.UTCDateTime(yr,mon,day,hr,mn,int(sec),int((sec-int(sec))*100000))
"""
##########################################################################
# get GG cat data
##########################################################################
cat = parse_ggcat('../../catalogue/data/GGcat-161025.csv')

##########################################################################
# loop thru events and get data
##########################################################################

for evnum, ev in enumerate(cat): 
    # only look for post 1990 data
    if ev['datetime'] >= datetime.datetime(1990, 1, 1, 0, 0):
        
        '''
        # for testing, get Moe data
        ev['datetime'] = datetime.datetime(2012,06,19,10,53)
        ev['lat'] = -38.304
        ev['lon'] = 146.200
        '''
        
        print evnum, ev['datetime']
        
        # convert datetime object to UTCdatetime
        dt = utcdatetime.UTCDateTime(ev['datetime'])
        #start_time=utcdatetime.UTCDateTime(yr,mon,day,hr,mn,int(sec),int((sec-int(sec))*100000))
       
        # Build event object
        evnt = Event(resource_id='GG_cat_' + str(evnum+1), creation_info='AU')
        
        origin = Origin()
        origin.time = ev['datetime']
        origin.longitude = ev['lon']
        origin.latitude = ev['lat']
        origin.depth = ev['dep']
        
        evnt.origins.append(origin)
        
        mag = Magnitude(creation_info='GG_cat')
        mag.mag = ev['prefmag']
        mag.magnitude_type = ev['prefmagtype']
        
        evnt.magnitudes.append(mag)
        

        ''' the time window to request the data will be 20 minutes, check maximum travel time and increase this value accordingly '''
        #end_time=start_time+960 # 16 minutes
        start_time = dt - datetime.timedelta(seconds=60)
        end_time   = dt + datetime.timedelta(seconds=960) # 16 minutes
        #end_time   = dt + datetime.timedelta(seconds=600) # 5 minutes
        
        
        ''' get all waveform data available, use wildcards to reduce the data volume and speedup the process,
        unfortunately we need to request few times for every number of characters that forms the station name '''
        # kluge to fix non-retrieval of data  - loop through alphabet integers
        for ch in range(ord('A'), ord('Z')+1):
            print 'Stations beginning with ', chr(ch)
            st_3 = client.get_waveforms("AU", chr(ch)+"??", "", "[BSEH]?[ENZ]", start_time,end_time)
            st_4 = client.get_waveforms("AU", chr(ch)+"???", "", "[BSEH]?[ENZ]", start_time,end_time)
            if ch == ord('A'):            
                if len(st_4) > 0:
                    st=st_3+st_4
                else:
                    st=st_3
            
            else:
                if len(st_4) > 0:
                    st+=st_3+st_4
                else:
                    st+=st_3

        # Cleanup duplicate traces returned by server
 #       st.merge(-1) #-1 method only merges overlapping or adjacent traces with same i            
        # Now sort the streams by station and channel
        st.sort()
        # Cleanup duplicate traces returned by server
        st.merge(1, fill_value=None) #1 method only merges overlapping or adjacent traces with same id
        # Now sort the streams by station and channel
        st.sort()
        
        # check if waves folder exists
        if not path.isdir('waves'):
            makedirs('waves')
            
        # set mseed filename
        msfile = path.join('waves', ev['datetime'].strftime('%Y%m%d%H%M')+'.mseed')
        
        # now write streams for each event to mseed
        st.write(msfile, format="MSEED")          
        
        forcecrash=blah
