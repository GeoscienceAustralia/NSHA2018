from obspy.core import UTCDateTime
from obspy.fdsn import Client
from obspy import read,Stream,Trace
import numpy as np
from obspy.core.util.geodetics import locations2degrees

def sn_test(tr,noise_win_len,sig_win_len,ratio):
    wl_10deg = 7*60. # see download_data_iris.py
    tp = tr.stats.starttime + wl_10deg
    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
    signal_tr = tr.copy().trim(starttime=tp-1,endtime=tp+sig_win_len)
    m1 = np.max(abs(noise_tr.data))
    m2 = np.max(abs(signal_tr.data))
    sn_test = 0 if (m1*ratio > m2) else 1
    return sn_test
def create_trace(st,sta_code,sta_coordinate,channel,eve_coordinate,eve_ot,freq4,ratio1,ratio2):
    tr = st.select(station=sta_code,channel = channel)
    dist = locations2degrees(eve_coordinate[0],eve_coordinate[1],
                                          sta_coordinate[0],sta_coordinate[1])
    if (len(tr)==1):
        tr = tr[0]
        tr.stats.data_available = 1
        if (hasattr(tr.stats,'response')):
            tr.stats.metadata_available = 1
            vel = tr.copy().detrend('demean').taper(0.05). \
                  remove_response(pre_filt=freq4,output="VEL")
            disp = tr.copy().detrend('demean').taper(0.05). \
                  remove_response(pre_filt=freq4,output="DISP")
            disp.resample(2.0)
            sig_win_len = 0.36 * dist * 111.194929703 + 60
            noise_win_len = 3 * 60
            tr.stats.sn_test1 = sn_test(vel,noise_win_len,sig_win_len,ratio1)
            tr.stats.sn_test2 = sn_test(disp,noise_win_len,sig_win_len,ratio2)
        else:
            tr.stats.metadata_available = 0
            tr.stats.sn_test1 = 0
            tr.stats.sn_test2 = 0
            
    else:
        tr = Trace()
        tr.stats.data_available = 0
        tr.stats.metadata_available = 0
        tr.stats.sn_test1 = 0
        tr.stats.sn_test2 = 0
        tr.trim(starttime=eve_ot-7*60,endtime=eve_ot+7*60,pad=True,fill_value=0.0) 
    
    tr.stats.distance = dist
    tr.stats.eve_coord = eve_coordinate
    tr.stats.eve_ot = eve_ot
    tr.stats["coordinates"] = {}
    tr.stats["coordinates"]["latitude"] = sta_coordinate[0]
    tr.stats["coordinates"]["longitude"] = sta_coordinate[1]
    return tr

e = np.array([9, -25.9665, 131.9755, 1.1, 2013, 6, 9, 14, 22, 12, 5.8, 5.43])
data_dir = './output/'
eve_id = 'eve_' + str(int(e[0]))
station_file = data_dir + 'stations_eve_' + str(int(e[0])) + '.txt'
waveform_file = data_dir + 'waveforms_eve_' + str(int(e[0])) + '.pkl'
eve_coordinate = (float(e[1]),float(e[2]))
eve_ot = UTCDateTime(int(e[4]),int(e[5]),int(e[6]),int(e[7]),int(e[8]),
                     float(e[9]))

st = read(waveform_file)
st_all = Stream()
freq4 = [0.01,0.03,15,20.]
ratio1 = 20.
ratio2 = 8.0
line = 'WC4 -19.9619 134.3397 6.38647400151 0 0 1 0 0 1'

sta = line.split(' ')
sta_code = sta[0]
sta_coordinate = (float(sta[1]),float(sta[2]))
tr_bhn = create_trace(st,sta_code,sta_coordinate,'BHN',eve_coordinate,
                      eve_ot,freq4,ratio1,ratio2)
tr_bhe = create_trace(st,sta_code,sta_coordinate,'BHE',eve_coordinate,
                      eve_ot,freq4,ratio1,ratio2)
tr_bhz = create_trace(st,sta_code,sta_coordinate,'BHZ',eve_coordinate,
                      eve_ot,freq4,ratio1,ratio2)
st_all.append(tr_bhn)
st_all.append(tr_bhe)
st_all.append(tr_bhz)

### test S/N tests
##bitrate = 2.0
##station = 'KMBL'
##channel = 'BHN'
##waveform = './output/waveforms_eve_1.pkl'
##distance = 3.67* 111.194929703
##sig_win_len = (0.36 * distance) + 60
##noise_win_len = 3*60
##ratio = 20.
##def sn_test1(tr,noise_win_len,sig_win_len,ratio):
##    wl_10deg = 7*60. # see download_data_iris.py
##    tp = tr.stats.starttime + wl_10deg
##    tr.detrend('demean')
##    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
##    signal_tr = tr.copy().trim(starttime=tp-1,endtime=tp+sig_win_len)
##    m1 = np.max(abs(noise_tr.data))
##    m2 = np.max(abs(signal_tr.data))
##    sn_test1 = 0 if (m1*ratio > m2) else 1
##    print m2/m1
##    print m1
##    print m2
##    return sn_test1
##
##st = read(waveform)
##tr = st.select(station=station,channel=channel)[0]
##vel = tr.copy()
##frq4 = [0.01,0.03,15,20.]
##vel.remove_response(pre_filt=frq4,output="VEL")
###disp1 = vel.copy().integrate()
##disp1 = tr.copy().remove_response(pre_filt=frq4,output="DISP")
##factor = int(disp1.stats.sampling_rate/bitrate)
##while (factor>16):
##    disp1.decimate(5)
##    factor = int(factor/5)
##disp1.decimate(factor)
##
###disp2 = vel.copy().integrate()
##disp2 = tr.copy().remove_response(pre_filt=frq4,output="DISP")
##disp2 = disp2.resample(bitrate)
##a = sn_test1(disp2,noise_win_len,sig_win_len,ratio)
##print a

### test downsampling
##bitrate = 2.0
##station = 'KMBL'
##channel = 'BHE'
##waveform = './output/waveforms_eve_1.pkl'
##st = read(waveform)
##tr = st.select(station=station,channel=channel)[0]
##vel = tr.copy()
##frq4 = [0.005,0.01,15,20.]
##vel.remove_response(pre_filt=frq4,output="VEL")
##disp1 = vel.copy().integrate()
##factor = int(disp1.stats.sampling_rate/bitrate)
##while (factor>16):
##    disp1.decimate(5)
##    factor = int(factor/5)
##disp1.decimate(factor)
##
##disp2 = vel.copy().resample(bitrate)
##client = Client('IRIS')
##
##
### test getting bulk data
##ot = UTCDateTime(2014,12,1,10,25,16)
##t1 = ot
##t2 = ot+7*60
##bulk = [('AU','AS31','*','BHZ',t1,t2),
##        ('AU','WC4','*','BHZ',t1,t2)]
##
##
##st = client.get_waveforms_bulk(bulk)


# plot
##import matplotlib.pyplot as plt
##from obspy import read
##import pdb
##
###st = read()
##st = read('test.pkl')
##st.normalize()
##ax = plt.axes()
##ylabels = []
##for i, tr in enumerate(st):
##    ylabels.append(tr.id)
##    data = tr.data
##    shift = i - data[0]
##    data2 = data + shift
##    ax.plot(tr.times(), data2, linestyle="-", color='black',
##        linewidth=1.5)
##    
##ax.set_yticks(range(len(st)))
##ax.set_yticklabels(ylabels)
##    
##plt.show()




