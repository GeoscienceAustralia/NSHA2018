
# coding: utf-8

# In[1]:

import local_magnitude
#from NSHA2018.magnitude.mw.Main.scripts import modules
from magnitude.mw.Main.scripts import modules


# In[2]:

#!/opt/anaconda-2.3.0/bin/python2
"""#!/opt/antelope/5.4/bin/python"""
#get_ipython().magic(u'matplotlib inline')
import numpy as np
from obspy.core import utcdatetime, event
from obspy.core.event import Catalog, Event, Magnitude, Origin, StationMagnitude
#from obspy.neic.client import Client
from obspy.clients.neic.client import Client
from obspy.io.xseed import Parser
#from obspy.core.util import gps2DistAzimuth
from obspy.geodetics import gps2dist_azimuth as gps2DistAzimuth
#from obspy.taup import TauPyModel
from obspy.signal import invsim as inv
from obspy.io.xseed.utils import SEEDParserException
# plot all traces into one pdf file
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot
from matplotlib.dates import date2num
# we will use dataless seed from IRIS to get station information
parser = Parser("../../data/AU.dataless")
# travel-time model will be iasp91 but could be any
from obspy import taup
vel_model = taup.TauPyModel(model="iasp91")
#from obspy.taup.TauPyModel import get_travel_times
# local modules
from magnitude.ml import local_magnitude
r_earth = 6371
def sind(x): return np.sin(x / 180. * np.pi)
def cosd(x): return np.cos(x / 180. * np.pi)
def tand(x): return np.tan(x / 180. * np.pi)
def arcsind(x): return np.arcsin(x) / np.pi * 180
def arccosd(x): return np.arccos(x) / np.pi * 180
def arctand(x): return np.arctan(x) / np.pi * 180
def gps2DistDegree(lat1, lon1, lat2, lon2):
    return arccosd(sind(lat1) * sind(lat2) +
                   cosd(lat1) * cosd(lat2) * cosd(lon1 - lon2))


# In[ ]:




# In[3]:

# here we define how we measure peak to peak amplitude and period
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


# In[4]:

# initialize the cwb port
client=Client(host='10.7.161.60',port=2061,debug=False)#, nonice=True)
eq=[]
# here we read all events line by line


# In[ ]:

# Instantiate catalogue object
catalogue = Catalog()
# Build Quakeml Event object


# In[ ]:

# file to store outputs
event_num = 0
f_out = open('amp_cat.txt', 'w')
with open("../../../GGCat_2014.csv",'r') as cat:
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
        lat=float(eq[13]) #Latitudes in old Gary Gibson catalogue given in southern hemisphere coordindates
        dep=float(eq[14])
        mag_type = eq[16]
        mag_pref = float(eq[17])
        mag_source = eq[31].split(' ')[-1]

        start_time=utcdatetime.UTCDateTime(yr,mon,day,hr,mn,int(sec),int((sec-int(sec))*100000))
#        ''' correct the time if not UTC '''
# Is this really needed??? Looks like may already be converted??
#        if code=="AEST" :
#            start_time-=36000
#        if code=="WITA" or code=="AWST":
#            start_time-=28800
        
        # Build event object
        event = Event(resource_id='GG_cat_' + str(event_num), creation_info='JG')
        event_num += 1
        origin = Origin()
        origin.time = start_time
        origin.longitude = lon
        origin.latitude = lat
        origin.depth = dep
        event.origins.append(origin)
        mag = Magnitude(creation_info='GG_cat')
        mag.mag = mag_pref
        mag.magnitude_type = mag_type
        event.magnitudes.append(mag)
        

        ''' the time window to request the data will be 20 minutes, check maximum travel time and increase this value accordingly '''
        end_time=start_time+1200 # 20 minutes - to catch ~Rayleigh waves at 20 degrees distance (assume v>=3km/s)
        ''' get all waveform data available, use wildcards to reduce the data volume and speedup the process,
        unfortunately we need to request few times for every number of characters that forms the station name '''
        st_3 = client.get_waveforms("AU", "???", "", "[BS]?[ENZ]", start_time,end_time)
        st_4 = client.get_waveforms("AU", "????", "", "[BS]?[ENZ]", start_time,end_time)
        if len(st_4) > 0:
            st=st_3+st_4
        else:
            st=st_3
        print st
        # Cleanup duplicate traces returned by server
 #       st.merge(-1) #-1 method only merges overlapping or adjacent traces with same i            
        # Now sort the streams by station and channel
        st.sort()
        # Cleanup duplicate traces returned by server
        st.merge(1, fill_value=None) #1 method only merges overlapping or adjacent traces with same id
        # Now sort the streams by station and channel
        st.sort()
        
        # Filter the stream - should it be done here on by trace???
        st.filter('bandpass', freqmin=0.5, freqmax=10., corners=6)#, df=sample_rate)
        # Demean
        st.detrend('demean')
  #      st.filter(bandpass, freqmin=0.5, freqmax=10., )
        # PDF file to plot all traces for this event
##        figpdf = PdfPages('test.pdf')
        try:
            st_figure = st.plot(equal_scale=False, handle=True)
        except:
            pass

        ''' first we will print the record of the earthquake as in catalogue '''
        print line
        
        counter = 0 # For debugging only
        mag_list = []
        distance_list = []
        for tr in st:
            counter +=1
            ''' get station ID'''
            seedid=tr.get_id()
            channel = tr.stats.channel
            # Filter to only get stations within 20 degrees
            try:
                ''' get station coordinates from dataless seed '''
                tr.stats.coordinates = parser.get_coordinates(seedid,start_time)
                tr.stats.distance = gps2DistDegree(tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,lat,lon)            
            except (SEEDParserException,AssertionError):
         #       print  tr.stats['station'],tr.stats['channel'],'-1,-1,-1'
                continue
            if tr.stats.distance > 20: # only use stations within 20 degrees
                continue
            else:
                pass
            

#            numobs = tr.count()
#            print 'count', numobs
#            print 'seedid', seedid
#            print 'sample rate', tr.stats.sampling_rate
            sample_rate = tr.stats.sampling_rate


            try:
                tr.stats.great_circle_distance, azf, azb =                     gps2DistAzimuth(tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,lat,lon)
                travel_times=vel_model.get_travel_times(dep, tr.stats.distance,                                                         phase_list = ['P', 'p', 'Pn', 'S', 'Sn'])#,model="iasp91") # need to make more efficient in order
                # to only calculate once per station
                #print travel_times, type(travel_times)
                #try:
                    #for arrival in travel_times:
                    #    print arrival.phase.name
                #   arrivals= (item for item in travel_times if item.phase.name=='S' or item.phase.name=='Sn').next()
                #except StopIteration:
                #   print "WARNING: reference station ",tr.stats['station']," does not have S or Sn phases"
                P_time = None
                S_time = None
                for arrival in travel_times:
                 #   print arrival.phase.name
                    if arrival.phase.name == 'P' or arrival.phase.name == 'p' or arrival.phase.name == 'Pn':
                        if P_time is not None:
                            P_time = min(P_time, arrival.time)
                        else:
                            P_time = arrival.time
                    if arrival.phase.name == 'S' or arrival.phase.name == 'Sn':
                        if S_time is not None:
                            S_time = min(S_time, arrival.time)
                        else:
                            S_time = arrival.time
#                print 'P_time', P_time
#                print 'S_time', S_time
                # Convert travel times to days and add start-time
                # See https://github.com/obspy/obspy/blob/master/obspy/imaging/waveform.py line 694
                # Then plot over waveforms
                SECONDS_PER_DAY = 3600.0 * 24.0
                try:
                    p_value = ((P_time / SECONDS_PER_DAY) +
                                date2num(tr.stats.starttime.datetime))
                except ValueError:
                    p_value = None
                try:
                    s_value = ((S_time / SECONDS_PER_DAY) +
                                date2num(tr.stats.starttime.datetime))
                except:
                    s_value = None
                    
                try:
                    st_figure.axes[counter-1].axvline(p_value, 0, 1, c='r')
                    st_figure.axes[counter-1].axvline(s_value, 0, 1, c='b')
                except:
                    pass
                # Save figures
             #   fig_name = 
               # pyplot.savefig()
                # Signal to noise tests
                try:
                    #sn = modules.sn_test(tr, start_time+P_time, 3*60., 0.36*tr.stats.great_circle_distance/1000. + 60, 20.)
                    sn, m1, m2 = modules.sn_test(tr, start_time+P_time, 30., 60, 5.)#3*(S_time-P_time), 5.)
                    snr = modules.snr_test(tr, start_time+P_time, 30., 60, [0.5, 10], 0.75)
              #  print 'sn', sn, m1, m2
              #  print snr
                except:
                    print 'Could not do signal to noise test'
                    sn = 0
                    snr = 0
               # snr = modules.snr_test()
#                print arrivals
#                for key, value in arrivals.iteritems():
#                        print key, value
##                S_time = start_time + arrivals['time']
##                figure = tr.plot(handle=True)
                # Add S arrival time
##                pyplot.scatter([S_time],[0], marker = 'o', s=40, c='r')
##                try:
##                   figpdf.savefig(figure)
##                    figure.close()
##               except:
##                    pass
               # figpdf.close()
#                  print travel_times
                ''' lets select 20 seconds window from S-wave p2p measurements '''
#               print start_time+arrivals['time'],start_time+arrivals['time']+20,end_time
 #               wave=tr.slice(start_time+arrivals['time'],start_time+arrivals['time']+20)
 #               if len(wave.data)>0:
 #                  amplitude,period,delay=max_p2t(wave.data,wave.stats.delta)
 #               else:
 #                   amplitude,period,delay = None,None,None
 #               if amplitude > 0. and period > 0. :
                paz=parser.get_paz(seedid,start_time)
                ''' calib is already applied therefore we set sensitivity to 1 '''
                paz['sensitivity']=1.
                # Define Wood Anderson paz based on Uhrhammer 1990 sensitivity of 2080
                # rather than theoretical 2800 as used by obpy
                # wa_amp=inv.estimateWoodAndersonAmplitude(paz,amplitude,period)
                paz_wa = {'sensitivity': 2080, 'zeros': [0j], 'gain': 1,
                          'poles': [-6.2832 - 4.7124j, -6.2832 + 4.7124j]}
                # Now convert to WA spectra
                wa_tr = tr.simulate(paz_remove=paz, paz_simulate=paz_wa)
                # Repeat signal to noise tests on displacement
                try:
                    sn_disp, m1, m2 = modules.sn_test(wa_tr, start_time+P_time, 30., 60., 5.)
                except:
                    sn_disp = 0
       #         print 'sn_dip', sn_disp, m1, m2
                # Plot WA record
##                figure = wa_tr.plot(handle=True)
##                try:
##                    figpdf.savefig(figure)
##                    figure.close()
##                except:
##                    pass
                # Calculate WA amplitudes and periods
#                wave=wa_tr.slice(start_time+arrivals['time'],start_time+arrivals['time']+20)
#                print len(wave.data)
                if len(wa_tr.data)>0:
                    wa_amp,period,delay=max_p2t(wa_tr.data,wa_tr.stats.delta)
                    local_mag = local_magnitude.calculate_local_magnitude(wa_amp/10e6,                                                                           [lon,lat,dep],                                                                          tr.stats.great_circle_distance/1000.)
                                        
                else:
                    wa_amp,period,delay = None,None,None
                    local_mag = None
                ''' now we print every station-component measurement '''
                #local_magnitude.calculate_local_magnitude
                
     #           print tr.stats['station'],tr.stats['channel'],tr.stats['calib'], \
      #              tr.stats.great_circle_distance/1000, wa_amp/10e6,period, local_mag
                mag_list.append(local_mag)
                distance_list.append(tr.stats.great_circle_distance/1000)
                outlist = [tr.stats['station'],tr.stats['channel'],tr.stats['calib'],                            tr.stats.coordinates.latitude,tr.stats.coordinates.longitude,                            lat, lon, dep, wa_amp, period, tr.stats.great_circle_distance,                            sn, snr, sn_disp]
                for item in outlist:
 #                   print item
                    f_out.write("%s," % item)
                f_out.write('\n')
                #else:
                #    print tr.stats['station'],tr.stats['channel'],tr.stats['calib'],amplitude,period
            except (SEEDParserException,AssertionError):
                pass
                #print  tr.stats['station'],tr.stats['channel'],'-1,-1,-1'
##            if counter > 10:
##                break
##    figpdf.close()
f_out.close()
#mean_mag = np.mean(mag_list)
##print 'mean ml', mean_mag
#pyplot.scatter(distance_list, mag_list)
#pyplot.plot([0,np.max(distance_list)], [mean_mag,mean_mag])
#pyplot.show()


# In[ ]:



