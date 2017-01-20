# Author: Hadi Ghasemi
# Date: 20-01-16

# Task: Download data within 1000km radius for events listed in catalogue.txt

# import modules
import pdb
import numpy as np
from obspy.fdsn import Client
from obspy.core import UTCDateTime
from obspy.core.util.geodetics import locations2degrees
from obspy.taup import TauPyModel

# default parameters
cat_file = './catalogue.txt'
client = Client('IRIS')
radius = 10. # in degrees
model = TauPyModel(model='iasp91')
wl_10deg = 7*60. # in second
output_dir = './output/'

# read catalogue.txt (see plan-phase1.txt for format details)
cat = np.loadtxt(cat_file)
cat = [cat[8]]
for line in cat:
    eve_lat = line[1]
    eve_lon = line[2]
    eve_depth = line[3]
    eve_ot = UTCDateTime(int(line[4]),int(line[5]),int(line[6]),int(line[7]),
                         int(line[8]),float(line[9]))
    print ("Downloading station inventory for eve_" + str(int(line[0])))
    inv = client.get_stations(network="AU",station="*",channel="BH*",starttime=eve_ot,
                              endtime=eve_ot+3600,latitude=eve_lat,longitude=eve_lon,
                              maxradius=radius)
    print ("completed for eve_" + str(int(line[0])))
    # bulk request from IRIS
    bulk_req = []
    S = []
    for sta in inv[0].stations:
        dist = locations2degrees(eve_lat,eve_lon,sta.latitude,sta.longitude)
        arrival = model.get_travel_times(source_depth_in_km=eve_depth,
                                       distance_in_degree=dist,phase_list=['P'])
        p_arrival = eve_ot + arrival[0].time
        t1 = p_arrival - wl_10deg
        t2 = p_arrival + wl_10deg
        bulk_req.append(('AU',sta.code,'*','BH*',t1,t2))
        s = sta.code + ' ' + str(sta.latitude) + ' ' + str(sta.longitude) \
            + ' ' + str(dist)
        S.append(s)
    print ("Downloading waveforms for eve_" + str(int(line[0])))
    st = client.get_waveforms_bulk(bulk_req, attach_response=True)
    print ("completed for eve_" + str(int(line[0])))
    # write stations.txt (see plan-phase1.txt for format details)
##    stations = []
##    for s in S:
##        sta_code = s.split(' ')[0]
##        tr_bhn = st.select(station=sta_code,channel ='BHN') 
##        tr_bhe = st.select(station=sta_code,channel ='BHE') 
##        tr_bhz = st.select(station=sta_code,channel ='BHZ') 
##        data_bhn = 1 if (len(tr_bhn)==1) else 0
##        data_bhe = 1 if (len(tr_bhe)==1) else 0
##        data_bhz = 1 if (len(tr_bhz)==1) else 0
##        if (data_bhn==1):
##            meta_data_bhn = 1 if (hasattr(tr_bhn[0].stats,'response')) else 0
##        else:
##            meta_data_bhn = 0
##        if (data_bhe==1):
##            meta_data_bhe = 1 if (hasattr(tr_bhe[0].stats,'response')) else 0
##        else:
##            meta_data_bhe = 0
##        if (data_bhz==1):
##            meta_data_bhz = 1 if (hasattr(tr_bhz[0].stats,'response')) else 0
##        else:
##            meta_data_bhz = 0
##        s = s + ' ' + str(data_bhn) + ' ' + str(data_bhe) + ' ' + str(data_bhz)\
##            + ' ' + str(meta_data_bhn) + ' ' + str(meta_data_bhe) + \
##            ' ' + str(meta_data_bhz)
##        stations.append(s)
    print (" writing waveforms into miniseed for eve_" + str(int(line[0])))
    output = output_dir + 'waveforms_eve_' + str(int(line[0])) + '.pkl'
    st.write(output, format = "PICKLE")
    print ("Finished!")
    print (" writing stations for eve_" + str(int(line[0])))
    output = output_dir + 'stations_eve_' + str(int(line[0])) + '.txt'
    f = open(output,'w')
    f.write("\n".join(map(lambda x: str(x), S)))
    f.close()
    print ("Finished!")
   
