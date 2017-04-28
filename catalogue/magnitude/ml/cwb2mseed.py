#!/usr/bin/env python
import numpy as np
import datetime
from os import path, makedirs
from obspy.core import utcdatetime, event
from obspy.core.event import Catalog, Event, Magnitude, Origin, StationMagnitude
import sys
from obspy.neic.client import Client
#from obspy.clients.neic.client import Client
from obspy.xseed import Parser
#from obspy.io.xseed import Parser
from obspy.xseed.utils import SEEDParserException
#from obspy.io.xseed.utils import SEEDParserException
from multiprocessing import Pool, Process, Queue, cpu_count
sys.path.append('../../')
from catalogue.parsers import parse_ggcat

# we will use dataless seed from IRIS to get station information
parser = Parser("../../data/AU.seed")


def sind(x): return np.sin(x / 180. * np.pi)
def cosd(x): return np.cos(x / 180. * np.pi)
def tand(x): return np.tan(x / 180. * np.pi)
def arcsind(x): return np.arcsin(x) / np.pi * 180
def arccosd(x): return np.arccos(x) / np.pi * 180
def arctand(x): return np.arctan(x) / np.pi * 180
def gps2DistDegree(lat1, lon1, lat2, lon2):
    return arccosd(sind(lat1) * sind(lat2) +
                   cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))


def get_sampling_rate(info,ev):

                channels = parser.getInventory()["channels"]
                time=utcdatetime.UTCDateTime(ev["datetime"])
                for channel in channels:
                    channel_id = channel.pop("channel_id")
                    net,sta,loc,comp=channel_id.split(".")
                    if not channel["end_date"]:
                        channel["end_date"]=utcdatetime.UTCDateTime.now()+1

                    if net == info["network"] and sta == info["station"] and loc == info["loc"] and comp==info["comp"]:
                        if utcdatetime.UTCDateTime(channel["start_date"]) <= time and utcdatetime.UTCDateTime(channel["end_date"]) >= time:
                           return int(channel["sampling_rate"])
                return int(0)


def get_cwb_data(ev):

     client=Client(host='10.7.161.60',port=2061,debug=False, timeout=60,nonice=True)

     origin = utcdatetime.UTCDateTime(ev['datetime'])


     ''' the time window to request the data will be 20 minutes, check maximum travel time and increase this value accordingly '''
     start_time = utcdatetime.UTCDateTime(origin - 60.)  # 1 minute
     end_time   = utcdatetime.UTCDateTime(origin + 960.) # 16 minutes

     ''' get all waveform data available, use wildcards to reduce the data volume and speedup the process,
     unfortunately we need to request few times for every number of characters that forms the station name '''

     st_3 = client.getWaveform("AU", "???", "", "[BS]??", start_time,end_time)
     st_4 = client.getWaveform("AU", "????", "", "[BS]??", start_time,end_time)

     if len(st_4) > 0:
        st=st_3+st_4
     else:
        st=st_3

     seed_sampling={}
     for tr in st:
            ''' get station ID'''
            seedid=str(tr.getId())


            ''' fix sampling rate if drift is detected '''

            info = {
                   "station":tr.stats.station,
                   "network":tr.stats.network,
                   "loc":tr.stats.location,
                   "comp":tr.stats.channel}

            ''' get sampling rate from metadata file '''
            sampling=get_sampling_rate(info,ev)
            if sampling:
                if sampling != tr.stats.sampling_rate :

                    if np.abs(sampling-tr.stats.sampling_rate) < 10:
                          tr.stats["sampling_rate"]=sampling
                    else:
                          print 'Station: ',seedid,' has sampling rate ',str(tr.stats.sampling_rate),' meanwhile ',str(sampling),' is assigned in dbmaster\n'
#           else:
#               print "Check your dataless seed completeness, no methadata for ",seedid," at ",start_time

            ''' fix sampling rate by interpolation if two streams have really big difference in sampling rate '''
            if seed_sampling.get(seedid,False):
                if seed_sampling[seedid]!=tr.stats.sampling_rate:
                   st_int=st.select(id=seedid)
                   max_samp=0
                   for trace in st_int:
                       if max_samp < trace.stats.sampling_rate:
                           max_samp = trace.stats.sampling_rate

                   for trace in st_int:
                       st.remove(trace)

                   st_int.resample(max_samp)

                   for trace in st_int:
                       st.append(trace)
                   print "Sampling mismatch found and corrected for ", seedid, " at ",start_time
            else:
                seed_sampling[seedid]=tr.stats.sampling_rate



     if not path.isdir('waves'):
            makedirs('waves')

     '''  set mseed filename '''

     msfile = path.join('waves', ev['datetime'].strftime('%Y%m%d%H%M')+'.mseed')

     st.sort()


     ''' check for masked arrays, there are some incompatibilities in obspy mseed libraries that should be fixed soon in feature versions. It is workaround for time being '''
     st.merge(method=0,fill_value='interpolate',interpolation_samples=0)


     for tr in st:
            if isinstance(tr.data,np.ma.masked_array):
                tr.data=tr.data.filled()
            tr.stats.mseed['encoding']='STEIM2'

     ''' now write streams for each event to mseed '''

     if len(st):
         st.write(msfile, format="MSEED")
         return int(1)
     else:
         return int(0)


def main(argv):


# Instantiate catalogue object
     catalogue = Catalog()

# Build Quakeml Event object

##########################################################################
# get GG cat data
##########################################################################

     cat = parse_ggcat('../../catalogue/data/test.csv')

##########################################################################
# loop thru events and get data
##########################################################################

     ''' Lets start to group the events by 10 and reminder '''
     istart=0
     chank=10
     for evnum, ev in enumerate(cat):
    # only look for post 1990 data
         if ev['datetime'] >= datetime.datetime(2010, 1, 1, 0, 0):
            istart=evnum
            break

     rem=(len(ev)-istart)%chank
     div=int((len(ev)-istart)/chank)
     pool_cpu=[]
     out_q = Queue()
     Abort_queue = False

     for i in xrange(div+1):
         i_start=i*chank
         i_stop=(i+1)*chank
         if i_stop > len(ev):
             i_stop=i_start+rem

         for j in xrange(i_start+istart,i_stop+istart):
             ev=cat[j-1]
             p = Process(target=get_cwb_data,args=(ev,))
             pool_cpu.append(p)
             p.start()

     for p in pool_cpu:
         if Abort_queue:
              p.terminate()
         p.join()


     if Abort_queue:
         exit(-1)




if __name__ == "__main__":
 main(sys.argv)

