# Author: Hadi Ghasemi
# Date: 27-01-16

# import modules
import pdb
import numpy as np
from obspy import read, Stream, Trace
from obspy.core import UTCDateTime
from obspy.core.util.geodetics import locations2degrees
import mtspec
from smtk.sm_utils import nextpow2


# default parameters
cat_file = './catalogue.txt'
data_dir = './output/'
freq4 = [0.01,0.03,15,20.]
ratio1 = 20.
ratio2 = 8.0
freq_range = [0.03,0.15]
perc = 0.75


# define modules
def snr_test(tr,noise_win_len,sig_win_len,freq_range,perc):
    wl_10deg = 7*60. # see download_data_iris.py
    tp = tr.stats.starttime + wl_10deg
    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
    signal_tr = tr.copy().trim(starttime=tp-1,endtime=tp+sig_win_len)
    noise_tr.detrend('demean').taper(0.05)
    signal_tr.detrend('demean').taper(0.05)
    n_fft = nextpow2(max(signal_tr.stats.npts,noise_tr.stats.npts))
    sig_spec,sig_freq = mtspec.mtspec(signal_tr.data,signal_tr.stats.delta,2,nfft=n_fft)
    sig_spec = np.sqrt(sig_spec)
    noise_spec,noise_freq = mtspec.mtspec(noise_tr.data,noise_tr.stats.delta,2,nfft=n_fft)
    noise_spec = np.sqrt(noise_spec)*(float(signal_tr.stats.npts)/noise_tr.stats.npts)
    sn_ratio = sig_spec/noise_spec
    index = np.where((sig_freq>=freq_range[0])&(sig_freq<=freq_range[1]))
    sn_ratio = sn_ratio[index]
    L1 = np.size(np.where(sn_ratio>=2.5))
    L2 = np.size(sn_ratio)
    snr_test = 1 if ((float(L1)/float(L2))>perc) else 0
    return snr_test
    
def sn_test(tr,noise_win_len,sig_win_len,ratio):
    wl_10deg = 7*60. # see download_data_iris.py
    tp = tr.stats.starttime + wl_10deg
    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
    signal_tr = tr.copy().trim(starttime=tp-1,endtime=tp+sig_win_len)
    m1 = np.max(abs(noise_tr.data))
    m2 = np.max(abs(signal_tr.data))
    sn_test = 0 if (m1*ratio > m2) else 1
    return sn_test
        
def create_trace(st,sta_code,sta_coordinate,channel,eve_coordinate,eve_ot,freq4,ratio1,ratio2,
                 freq_range,perc):
    tr = st.copy().select(station=sta_code,channel = channel)
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
            pdb.set_trace()
            sig_win_len = 0.36 * dist * 111.194929703 + 60
            noise_win_len = 3 * 60
            tr.stats.sn_test1 = sn_test(vel,noise_win_len,sig_win_len,ratio1)
            tr.stats.sn_test2 = sn_test(disp,noise_win_len,sig_win_len,ratio2)
            tr.stats.snr_test = snr_test(vel,noise_win_len,sig_win_len,freq_range,perc) 
        else:
            tr.stats.metadata_available = 0
            tr.stats.sn_test1 = 0
            tr.stats.sn_test2 = 0
            tr.stats.snr_test = 0 
            
    else:
        tr = Trace()
        tr.stats.data_available = 0
        tr.stats.metadata_available = 0
        tr.stats.sn_test1 = 0
        tr.stats.sn_test2 = 0
        tr.stats.snr_test = 0 
        tr.trim(starttime=eve_ot-7*60,endtime=eve_ot+7*60,pad=True,fill_value=0.00001)
        tr.stats.station = sta_code
        tr.stats.channel = channel
    
    tr.stats.distance = dist
    tr.stats.eve_coord = eve_coordinate
    tr.stats.eve_ot = eve_ot
    tr.stats["coordinates"] = {}
    tr.stats["coordinates"]["latitude"] = sta_coordinate[0]
    tr.stats["coordinates"]["longitude"] = sta_coordinate[1]
    return tr

#
cat = np.loadtxt(cat_file)

for e in cat:
    eve_id = 'eve_' + str(int(e[0]))
    station_file = data_dir + 'stations_eve_' + str(int(e[0])) + '.txt'
    waveform_file = data_dir + 'waveforms_eve_' + str(int(e[0])) + '.pkl'
    eve_coordinate = (float(e[1]),float(e[2]))
    eve_ot = UTCDateTime(int(e[4]),int(e[5]),int(e[6]),int(e[7]),int(e[8]),
                         float(e[9]))
    with open(station_file) as f:
        lines = f.readlines()
    st = read(waveform_file)
    st_all = Stream()   
    for line in lines:
        sta = line.split(' ')
        sta_code = sta[0]
        sta_coordinate = (float(sta[1]),float(sta[2]))
        tr_bhn = create_trace(st,sta_code,sta_coordinate,'BHN',eve_coordinate,
                              eve_ot,freq4,ratio1,ratio2,freq_range,perc)
        tr_bhe = create_trace(st,sta_code,sta_coordinate,'BHE',eve_coordinate,
                              eve_ot,freq4,ratio1,ratio2,freq_range,perc)
        tr_bhz = create_trace(st,sta_code,sta_coordinate,'BHZ',eve_coordinate,
                              eve_ot,freq4,ratio1,ratio2,freq_range,perc)
        st_all.append(tr_bhn)
        st_all.append(tr_bhe)
        st_all.append(tr_bhz)
    
    print (" writing waveforms into miniseed for " + eve_id)
    output = data_dir + 'processed_waveforms_' + eve_id + '.pkl'
    st_all.write(output, format = "PICKLE")
    print ("Finished!")

     

    
  
