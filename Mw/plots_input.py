# Author: Hadi Ghasemi
# Date: 21-01-16

# import modules
import pdb
import numpy as np
from obspy import read, Stream
from obspy.core import UTCDateTime
import os
import matplotlib.pyplot as plt


#
cat_file = './catalogue.txt'
data_dir = './output/'
output_plot = './output/plots/'


#
cat = np.loadtxt(cat_file)
for e in cat:
    eve_id = 'eve_' + str(int(e[0]))
    waveform_file = data_dir + 'processed_waveforms_' + eve_id + '.pkl'
    st_all = read(waveform_file)
    st_all.detrend('demean')
    st_all.sort(['distance'])
    st_all.normalize()
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()
    ylabels = []
    for i, tr in enumerate(st_all):
        tr.data[np.isnan((tr.data))]=0 # normalize traces without data!
        ylabels.append(tr.stats.station + '.' + tr.stats.channel)
        shift = i - tr.data[0]
        data = tr.data + shift
        times = [(tr.stats.starttime + t).datetime for t in tr.times()]
        color = 'yellow'
##        if (tr.stats.data_available==0)|(tr.stats.metadata_available==0):
##            color = 'black'
##        elif (tr.stats.sn_test1==0)&(tr.stats.sn_test2==0):
##            color = 'red'
##        elif (tr.stats.sn_test1==1)&(tr.stats.sn_test2==0):
##            color = 'blue'
##        elif (tr.stats.sn_test1==1)&(tr.stats.sn_test2==1):
##            color = 'green'
        if (tr.stats.snr_test==1):
            color = 'black'
        
        ax.plot(times, data, linestyle="-", color=color,
            linewidth=1.5,alpha=0.5)
    ax.set_yticks(range(len(st_all)))
    ax.set_yticklabels(ylabels,fontsize=10)
    print ("save the plot for " + eve_id)
    plt.savefig(output_plot + eve_id + '/' + eve_id + '_all_processed_data_snr.png')
   

##cat = np.loadtxt(cat_file)
##for e in cat:
##    eve_id = 'eve_' + str(int(e[0]))
##    station_file = data_dir + 'stations_eve_' + str(int(e[0])) + '.txt'
##    waveform_file = data_dir + 'waveforms_eve_' + str(int(e[0])) + '.pkl'
##    eve_coordinate = (float(e[1]),float(e[2]))
##    eve_ot = UTCDateTime(int(e[4]),int(e[5]),int(e[6]),int(e[7]),int(e[8]),
##                         float(e[9]))
##    with open(station_file) as f:
##        lines = f.readlines()
##    st = read(waveform_file)
##    st_all = Stream()
##    st_bhn = Stream()
##    st_bhe = Stream()
##    st_bhz = Stream()
##    for line in lines:
##        sta = line.split(' ')
##        sta_code = sta[0]
##        sta_lat = float(sta[1])
##        sta_lon = float(sta[2])
##        sta_dist = float(sta[3])* 111.194929703 * 1000.
##        data_bhn = int(sta[4])
##        data_bhe = int(sta[5])
##        data_bhz = int(sta[6])
##        if (data_bhn==1):
##            tr_bhn = st.copy().select(station=sta_code,channel ='BHN')[0]
##            tr_bhn.stats.distance = sta_dist
##            tr_bhn.stats["coordinates"] = {}
##            tr_bhn.stats["coordinates"]["latitude"] = sta_lat
##            tr_bhn.stats["coordinates"]["longitude"] = sta_lon
##            tr_bhn.trim(starttime=eve_ot,endtime=eve_ot+600)
##            st_bhn.append(tr_bhn)
##            st_all.append(tr_bhn)
##
##        if (data_bhe==1):
##            tr_bhe = st.copy().select(station=sta_code,channel ='BHE')[0]
##            tr_bhe.stats.distance = sta_dist
##            tr_bhe.stats["coordinates"] = {}
##            tr_bhe.stats["coordinates"]["latitude"] = sta_lat
##            tr_bhe.stats["coordinates"]["longitude"] = sta_lon
##            tr_bhe.trim(starttime=eve_ot,endtime=eve_ot+600)
##            st_bhe.append(tr_bhe)
##            st_all.append(tr_bhe)
##    
##        if (data_bhz==1):
##            tr_bhz = st.copy().select(station=sta_code,channel ='BHZ')[0]
##            tr_bhz.stats.distance = sta_dist
##            tr_bhz.stats["coordinates"] = {}
##            tr_bhz.stats["coordinates"]["latitude"] = sta_lat
##            tr_bhz.stats["coordinates"]["longitude"] = sta_lon
##            tr_bhz.trim(starttime=eve_ot,endtime=eve_ot+600)
##            st_bhz.append(tr_bhz)
##            st_all.append(tr_bhz)
##    if not os.path.exists(output_plot + eve_id):
##        os.makedirs(output_plot + eve_id)
##    st_bhn.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhn_km.png',
##                type='section',plot_dx = 100000)       
##        
##    st_bhn.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhn_deg.png',
##                type='section',dist_degree=True,ev_coord=eve_coordinate,
##                plot_dx = 1)
##   
##    st_bhe.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhe_km.png',
##                type='section',plot_dx = 100000)       
##        
##    st_bhe.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhe_deg.png',
##                type='section',dist_degree=True,ev_coord=eve_coordinate,
##                plot_dx = 1)
##    st_bhz.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhz_km.png',
##                type='section',plot_dx = 100000)       
##        
##    st_bhz.plot(outfile = output_plot + eve_id + '/' + eve_id + '_bhz_deg.png',
##                type='section',dist_degree=True,ev_coord=eve_coordinate,
##                plot_dx = 1)
##
##    st_all.sort(['distance'])
##    st_all.normalize()
##    fig = plt.figure(figsize=(10,10))
##    ax = plt.axes()
##    ylabels = []
##    for i, tr in enumerate(st_all):
##        ylabels.append(tr.id[3::])
##        shift = i - tr.data[0]
##        data = tr.data + shift
##        times = [(tr.stats.starttime + t).datetime for t in tr.times()]
##        ax.plot(times, data, linestyle="-", color='black',
##            linewidth=1.5,alpha=0.5)
##    ax.set_yticks(range(len(st_all)))
##    ax.set_yticklabels(ylabels,fontsize=10)
##    plt.savefig(output_plot + eve_id + '/' + eve_id + '_all_data.png')    
##    
