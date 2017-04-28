import pdb
import os
import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt
from obspy import read, read_inventory, Stream
from obspy.core import UTCDateTime
from obspy.core.event import Catalog, Event, Origin, Magnitude
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth, locations2degrees
from obspy.taup import TauPyModel
import shapefile
from shapely.geometry import Polygon, Point



##############################catalogue###########################################
def read_cat_ref(cat_file):
    """
    Parses a given refrence catalogue (in ascii format,see the header for details)
    output is Obspy catalogue object
    """
    cat_ref = np.loadtxt(cat_file,delimiter=',',skiprows=1)
    cat = Catalog()
    for i,e in enumerate(cat_ref):
        event = Event(resource_id='smi:local/='+str(i),creation_info='HG')
        origin = Origin()
        origin.time = UTCDateTime(int(e[2]),int(e[3]),int(e[4]),
                                  int(e[7]),int(e[8]),e[9])
        origin.longitude = e[0]
        origin.latitude = e[1]
        origin.depth = e[6] * 1000. #in meters
        event.origins.append(origin)
        if ~(np.isnan(e[10])):
            mag = Magnitude(creation_info='HER')
            mag.mag = e[10]
            mag.magnitude_type = 'Mw'
            event.magnitudes.append(mag)
        if ~(np.isnan(e[11])):
            mag = Magnitude(creation_info='MAR')
            mag.mag = e[11]
            mag.magnitude_type = 'Mw' 
            event.magnitudes.append(mag)
        if ~(np.isnan(e[12])):
            mag = Magnitude(creation_info='SIP')
            mag.mag = e[12]
            mag.magnitude_type = 'Mw' 
            event.magnitudes.append(mag)
        cat.append(event)
    return cat

def get_cat_iris(sc):
    """
    Gets the event catalogue from IRIS based on sc:search_criteria
    output is Obspy catalogue object
    """
    client = Client("IRIS")
    cat = client.get_events(starttime=sc['start_time'],endtime=sc['end_time'],
                            minmagnitude=sc['min_mag'],
                            minlatitude=sc['min_lat'], maxlatitude=sc['max_lat'],
                            minlongitude=sc['min_lon'], maxlongitude=sc['max_lon'],
                            orderby=sc['order'])
    return cat

def get_au_poly(shp_file):
    """
    extracts the polygon of the mainland from shp file
    output is Shapely polygon object
    """
    r = shapefile.Reader(shp_file)
    shapes = r.shapes()
    records = r.records()
    for record, shape in zip(records,shapes):
        lons,lats = zip(*shape.points)
        data = np.array((lons, lats)).T
        segs = []
        for i in range(1,len(shape.parts)):
            index = shape.parts[i-1]
            index2 = shape.parts[i]
            segs.append(data[index:index2])
        segs.append(data[index2:])
        A = np.zeros((len(segs),2))
        k = -1
        for seg in segs:
            k = k+1
            A[k] = [k,len(seg)]
    idx = np.where(A[:,1]==np.max(A[:,1]))[0]
    x = segs[idx][:,0]
    y = segs[idx][:,1]
    poly = Polygon(zip(x,y))
    return poly

def filter_cat_poly(cat,poly):
    """
    Gets the events within a polygon
    output is Obspy catalogue object
    """
    cat_filt = Catalog()
    for eve in cat:
        if poly.contains(Point(eve.origins[0].longitude,
                                eve.origins[0].latitude)):
            cat_filt.append(eve)
    return cat_filt

def merge_cats(cat1,cat2,dt,dl):
    """
    find events in cat1 which are also in cat2, then adds the magnitudes from cat1
    into cat2
    output is Obspy catalogue object
    """
    for eve in cat1:
        tc1 = 'time > ' + str(eve.origins[0].time-5.0)
        tc2 = 'time < ' + str(eve.origins[0].time+5.0)
        c = cat2.filter(tc1,tc2)
        cat_com = cat2.copy()
        for e in c:
            dist = gps2dist_azimuth(eve.origins[0].latitude,
                                     eve.origins[0].longitude,
                                     e.origins[0].latitude,
                                     e.origins[0].longitude)[0]
            if dist<10:
                cat_com.events.remove(e)
                for m in eve.magnitudes:
                    e.magnitudes.append(m)
                cat_com.append(e)
    return cat_com     

##############################preparation###########################################                
def get_sta_wf_iris(sp,cat):
    """
    Gets the station inventory and waveforms for events in cat.
    output writes station inventory as xml files for each event
    output writes waveforms as pkl files for each event
    """
    client = Client("IRIS")
    for num,eve in enumerate(cat):
        eve_lat = eve.origins[0].latitude
        eve_lon = eve.origins[0].longitude
        eve_ot = eve.origins[0].time
        eve_id = eve.resource_id.id.split('=')[1] 
        print (str(num) + " Downloading station inventory for eve_" + eve_id)
        inv = client.get_stations(network=sp['network'],station=sp['station'],
                                  channel=sp['channel'],starttime=eve_ot,
                                  endtime=eve_ot+3600,latitude=eve_lat,longitude=eve_lon,
                                  maxradius=sp['maxradius'])
        print ("completed for eve_" + eve_id)
        
        # bulk request from IRIS
        bulk_req = []
        for i,net in enumerate(inv):
            for sta in inv[i].stations:
                t1 = eve_ot - sp['t_before_ot']
                t2 = eve_ot + sp['t_after_ot']
                bulk_req.append((net.code,sta.code,'*',sp['channel'],t1,t2))
        print (str(num) + " Downloading waveforms for eve_" + eve_id)
        st = client.get_waveforms_bulk(bulk_req, attach_response=True)
        print ("completed for eve_" + eve_id)
    
        print ("writing waveforms into pkl file for eve_" + eve_id)
        output = sp['output_dir_wf'] + 'waveforms_eve_' + eve_id + '.pkl'
        st.write(output, format = "PICKLE")
        print ("completed for eve_" + eve_id)
        print ("writing stations for eve_" + eve_id)
        output = sp['output_dir_sta'] + 'stations_eve_' + eve_id + '.xml'
        inv.write(output,format='STATIONXML')
        print ("completed for eve_" + eve_id)


def todisplacement(tr,freq4,bitrate):
    disp = tr.copy().detrend('demean').taper(0.05). \
           remove_response(pre_filt=freq4,output="DISP")
    disp.resample(bitrate)
##    # mouse-trp algorithm
##    factor = int(disp.stats.sampling_rate / bitrate)
##    while (factor > 16):
##        disp.decimate(5)
##        factor = int(factor / 5)
##    disp.decimate(factor)
    return disp

def est_p_time(depth,dist,ot):

        model = TauPyModel(model='iasp91')
        arrival = model.get_travel_times(source_depth_in_km=depth,
                                           distance_in_degree=dist,phase_list=['p','P'])
        p_arrival = ot + arrival[0].time
        return p_arrival

def tr_fft(tr):
    T = tr.stats.delta
    N = tr.stats.npts
    s = scipy.fftpack.fft(tr.data)
    s = 2.0/N * np.abs(s[:N/2])
    f = np.linspace(0.0, 1.0/(2.0*T), N/2)
    return s,f
    

   
def sn_test(tr,tp,noise_win_len,sig_win_len,ratio):
    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
    signal_tr = tr.copy().trim(starttime=tp-2,endtime=tp+sig_win_len)
    m1 = np.max(abs(noise_tr.data))
    m2 = np.max(abs(signal_tr.data))
    sn_test = 0 if (m1*ratio > m2) else 1
    return sn_test, m1, m2



def snr_test(tr,tp,noise_win_len,sig_win_len,freq_range,perc):
    noise_tr = tr.copy().trim(starttime=tp-noise_win_len-5,endtime=tp-5)
    signal_tr = tr.copy().trim(starttime=tp-1,endtime=tp+sig_win_len+1)
    noise_tr.detrend('demean').taper(0.05)
    signal_tr.detrend('demean').taper(0.05)
    sig_spec,sig_freq = tr_fft(signal_tr)
    noise_spec,noise_freq = tr_fft(noise_tr)
    freq = np.linspace(freq_range[0],freq_range[1],500)
    sig_spec2 = np.interp(freq,sig_freq,sig_spec)
    noise_spec2 = np.interp(freq,noise_freq,noise_spec)
    sn_ratio = sig_spec2/noise_spec2
    L1 = np.size(np.where(sn_ratio>=2.5))
    L2 = np.size(sn_ratio)
    snr_test = 1 if ((float(L1)/float(L2))>perc) else 0
    return snr_test

def process_wf_tr(tr,eve,sta,proc_params):
    """
    Gets the trace (1 channel, raw data) for 1 (eve,sta) pair
    output trace with added processing results to the header 
    
    """
    freq4 = proc_params['freq4']
    bitrate = proc_params['bitrate']
    ratio_vel = proc_params['ratio_vel']
    ratio_disp = proc_params['ratio_disp']
    freq_range_snr = proc_params['freq_range_snr']
    per_snr = proc_params['per_snr']
    vel = tr.copy().detrend('demean').taper(0.05). \
                  remove_response(pre_filt=freq4,output="VEL")
    disp = todisplacement(tr,freq4,bitrate)
    dist = locations2degrees(eve.origins[0].latitude,
                             eve.origins[0].longitude,
                             sta.latitude,sta.longitude)
    tp = est_p_time(eve.origins[0].depth/1000.,
                    dist,eve.origins[0].time)
    sig_win_len = 0.36 * dist * 111.194929703 + 60
    sn_test1 = sn_test(vel,tp,3*60.,sig_win_len,ratio_vel)
    sn_test2 = sn_test(disp,tp,3*60.,sig_win_len,ratio_disp)
    sn_test3 = snr_test(vel,tp,sig_win_len,
                        sig_win_len,freq_range_snr,per_snr)
    
    tr_proc = tr.copy()
    tr_proc.stats.eve_coord = (eve.origins[0].latitude,eve.origins[0].longitude)
    tr_proc.stats.eve_ot = eve.origins[0].time
    tr_proc.stats["coordinates"] = {}
    tr_proc.stats["coordinates"]["latitude"] = sta.latitude
    tr_proc.stats["coordinates"]["longitude"] = sta.longitude
    tr_proc.stats.distance = dist 
    tr_proc.stats.p_time = tp   
    tr_proc.stats.sn_test1 = sn_test1
    tr_proc.stats.sn_test2 = sn_test2
    tr_proc.stats.sn_test3 = sn_test3

    return tr_proc

    
def process_wf_eve(eve,inv,st,proc_params):
    st_proc = Stream()
    for i,net in enumerate(inv):
        for sta in inv[i].stations:
            st_sta = st.copy().select(network=net.code, station=sta.code)
            st_sta.merge(method=0,fill_value='latest')
            for tr in st_sta:
                tr_proc = process_wf_tr(tr,eve,sta,proc_params)
                st_proc.append(tr_proc)
    return st_proc

def plot_proc_wf(st,eve_id,plot_dir):
    st.detrend('demean')
    st.sort(['distance'])
    st.normalize()
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes()
    ylabels = []
    for i, tr in enumerate(st):
        ylabels.append(tr.stats.network + '.' + \
                       tr.stats.station + '.' + tr.stats.channel)
        shift = i - tr.data[0]
        data = tr.data + shift
        times = [(tr.stats.starttime + t).datetime for t in tr.times()]
        color = 'red'
        if (tr.stats.sn_test1==1)|(tr.stats.sn_test2==1)|\
           (tr.stats.sn_test3==1):
            color = 'green'
        ax.plot(times, data, linestyle="-", color=color,
            linewidth=1.5,alpha=0.5)
    ax.set_yticks(range(len(st)))
    ax.set_yticklabels(ylabels,fontsize=10)
    plt.savefig(plot_dir + eve_id + '_processed_wf.png')

def process_wf_cat(cat,sta_dir,wf_dir,proc_params,plot_dir):
    for eve in cat:
        eve_id = eve.resource_id.id.split('=')[1]
        inv_file = sta_dir+'stations_eve_' + eve_id + '.xml'
        wf_file = wf_dir+'waveforms_eve_' + eve_id + '.pkl'
        inv = read_inventory(inv_file)
        st = read(wf_file)
        print ("processing eve_" + eve_id)
        st_proc = process_wf_eve(eve,inv,st,proc_params)
        print ("writing processed waveforms into pkl file for eve_" + eve_id)
        st_proc.write(wf_file, format = "PICKLE")
        print ("completed for eve_" + eve_id)
        print ("plotting processed waveforms into png file for eve_" + eve_id)
        plot_proc_wf(st_proc,eve_id,plot_dir)
        print ("completed for eve_" + eve_id)        

    

def write_kiwi_input_file(e,proc_params,kiwi_dir,prefix):
    # general params
    eve_id = prefix + e.resource_id.id.split('=')[1]
    inp_file = kiwi_dir['kiwi_work_dir'] + 'rapidinv.inp.' + eve_id
    data_dir_eve = kiwi_dir['kiwi_data_dir'] + eve_id
    result_dir_eve = kiwi_dir['kiwi_result_dir'] + eve_id
    gfdb_dir = kiwi_dir['kiwi_gfdb_dir']
    confidence_int = 95
    epi_dist_min = 0
    epi_dist_max = 1000
    num_bootstrap = 200
    num_inv_steps = 1
    station_inp_file = 'stations.dat'
    sw_weight_dist = 'False'
    sw_filternoisy = 'True'
    # inversion step 1
    f_c = proc_params['freq_range_snr']
    depth_upperlim = 0.0
    depth_bottomlim = 30.
    inv_mode_step1 = 'invert_dmsdsok'
    loops_sds_conf = 1
    phases_to_use_st1 = 'a'
    weight_a_st1 = 1
    win_start_a_st1 = 0.02
    # writing
    data = []
    data.append('# GENERAL PARAMETERS' + '\n')
    data.append('GFDB_STEP1       ' + gfdb_dir + '\n')
    data.append('INVERSION_DIR    ' + result_dir_eve + '\n')
    data.append('DATA_DIR         ' + data_dir_eve + '\n')
    data.append('DATA_FILE        ' + eve_id + '\n')
    data.append('LATITUDE_NORTH   ' + \
                str(e.origins[0].latitude) + '\n')
    data.append('LONGITUDE_EAST   ' + \
                str(e.origins[0].longitude) + '\n')
    data.append('YEAR             ' + \
                str(e.origins[0].time.year) + '\n')
    data.append('MONTH            ' + \
                str(e.origins[0].time.month)+ '\n')
    data.append('DAY              ' + \
                str(e.origins[0].time.day)+ '\n')
    data.append('HOUR             ' + \
                str(e.origins[0].time.hour)+ '\n')
    data.append('MIN              ' + \
                str(e.origins[0].time.minute)+ '\n')
    data.append('SEC              ' + \
                str(e.origins[0].time.second)+ '\n')
    data.append('CONFIDENCE_INT   ' + str(confidence_int) + '\n')
    data.append('EPIC_DIST_MIN    ' + str(epi_dist_min) + '\n')
    data.append('EPIC_DIST_MAX    ' + str(epi_dist_max) + '\n')
    data.append('NUM_BOOTSTRAP    ' + str(num_bootstrap) + '\n')
    data.append('NUM_INV_STEPS    ' + str(num_inv_steps) + '\n')
    data.append('STAT_INP_FILE    ' + station_inp_file + '\n')
    data.append('SW_WEIGHT_DIST   ' + sw_weight_dist + '\n')
    data.append('SW_FILTERNOISY   ' + sw_filternoisy + '\n')
    
    data.append('\n' + '# INVERSION STEP 1 PARAMETERS' + '\n')
    data.append('BP_F1_STEP1      ' + str(f_c[0]) + '\n')
    data.append('BP_F2_STEP1      ' + str(f_c[0]) + '\n')
    data.append('BP_F3_STEP1      ' + str(f_c[1]) + '\n')
    data.append('BP_F4_STEP1      ' + str(f_c[1]) + '\n')
    data.append('DEPTH_UPPERLIM   ' + str(depth_upperlim) + '\n')
    data.append('DEPTH_BOTTOMLIM   ' + str(depth_bottomlim) + '\n')
    data.append('INV_MODE_STEP1    ' + inv_mode_step1 + '\n')
    data.append('LOOPS_SDS_CONF    ' + str(loops_sds_conf) + '\n')
    data.append('PHASES_2_USE_ST1 ' + phases_to_use_st1 + '\n')
    data.append('WEIGHT_A_ST1     ' + str(weight_a_st1) + '\n')
    data.append('WIN_START_A_ST1  ' + str(win_start_a_st1) + '\n')

    with open(inp_file,'w') as file:
        data = file.writelines(data)
        
def write_kiwi_station_file(inv,kiwi_dir,eve_id):
    output_dir = kiwi_dir['kiwi_data_dir']+ eve_id
    if not os.path.exists(output_dir):
        os.makedirs(output_dir) 
    output_file = output_dir + '/stations.dat'
    data = []
    k = -1
    for i,net in enumerate(inv):
        for sta in inv[i].stations:
            k = k + 1
            data.append(str(int(k)) + ' ' +  sta.code + ' ' + \
                        str(sta.latitude) + ' ' + str(sta.longitude) + '\n')       
    with open(output_file,'w') as file:
        data = file.writelines(data)
        
def write_kiwi_disp_wf(st,proc_params,kiwi_dir,eve_id):
    freq4 = proc_params['freq4']
    bitrate = proc_params['bitrate']
    for tr in st:
        if (tr.stats.sn_test1==1)|(tr.stats.sn_test2==1)|\
           (tr.stats.sn_test3==1):
##        if (tr.stats.network=='AU'):
            disp = todisplacement(tr,freq4,bitrate)
            disp_file = kiwi_dir['kiwi_data_dir'] + eve_id + \
                        '/DISPL.'+ tr.stats.station + '.' + \
                        tr.stats.channel
            disp.write(disp_file,format='MSEED')
            

    
def run_inv_eve(eve,wf,inv,proc_params,kiwi_dir,output_pre):
     current_dir = os.getcwd()
     eve_id = output_pre + eve.resource_id.id.split('=')[1]
     write_kiwi_input_file(eve,proc_params,kiwi_dir,output_pre)
     write_kiwi_station_file(inv,kiwi_dir,eve_id)
     write_kiwi_disp_wf(wf,proc_params,kiwi_dir,eve_id)
     os.chdir(kiwi_dir['kiwi_work_dir'])
     os.system('python rapidinv12.py' + ' ' + 'rapidinv.inp.' + eve_id)
     os.chdir(current_dir)

def parse_kiwi_output(kiwi_dir,eve_id,output_pre):
    kiwi_output_file = kiwi_dir['kiwi_result_dir'] + output_pre + eve_id \
                       + '/step1.earthquakeinfo.dat'
    if os.path.isfile(kiwi_output_file):
        with open(kiwi_output_file) as f:
            lines = f.readlines()
        Mw = float(lines[23].split(' ')[-1])
        depth = float(lines[25].split(' ')[-1][:-3])
        mis_fit = float(lines[29].split(' ')[-1])
        no_traces = int(lines[39].split(' ')[6])
    else:
        Mw,depth,mis_fit = np.array([np.nan,np.nan,np.nan])
    return np.array([Mw,mis_fit,depth])

def run_inv_cat(cat,wf_dir,inv_dir,proc_params,kiwi_dir,output_pre,results_dir):
    results = np.zeros((len(cat),3))
    for i,eve in enumerate(cat):
        eve_id = eve.resource_id.id.split('=')[1]
        inv_file = inv_dir+'stations_eve_' + eve_id + '.xml'
        wf_file = wf_dir+'waveforms_eve_' + eve_id + '.pkl'
        wf = read(wf_file)
        inv = read_inventory(inv_file)
        run_inv_eve(eve,wf,inv,proc_params,kiwi_dir,output_pre)
        results[i,:] = parse_kiwi_output(kiwi_dir,eve_id,output_pre)
    result_file = results_dir + 'results_inv.txt'
    header = 'Mw_inv misfit depth'
    np.savetxt(result_file,results,fmt='%.2f',header=header)
    
        
    
        

