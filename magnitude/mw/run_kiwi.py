# Author: Hadi Ghasemi
# Date: 12-02-16

# Task: run kiwi tools

# import modules
import pdb
import numpy as np
from obspy import read
import os
import matplotlib.pyplot as plt
from operator import sub
import matplotlib.patches as mpatch
from matplotlib.collections import PatchCollection

# default parameters
cat_file = './catalogue.txt'
data_dir = './output/'
kiwi_work_dir = '../WORK/'
kiwi_data_dir = '../DATA/'
kiwi_result_dir = '../RESULTS/'
current_dir = os.getcwd()

freq4 = [0.01,0.03,15,20.]
bitrate = 2.0

#  define modules
def write_kiwi_input_file(e,kiwi_work_dir):
    # general params
    eve_id = 'auto_eve_' + str(int(e[0]))
    inp_file = kiwi_work_dir + 'rapidinv.inp.' + eve_id
    confidence_int = 95
    epi_dist_min = 0
    epi_dist_max = 1000
    num_bootstrap = 200
    num_inv_steps = 1
    station_inp_file = 'stations.dat'
    sw_weight_dist = 'False'
    sw_filternoisy = 'True'
    # inversion step 1
    f_c = [0.035,0.035,0.15,0.15]
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
    data.append('INVERSION_DIR    ../RESULTS/' + eve_id + '\n')
    data.append('DATA_DIR         ../DATA/' + eve_id + '\n')
    data.append('DATA_FILE        ' + eve_id + '\n')
    data.append('LATITUDE_NORTH   ' + str(float(e[1])) + '\n')
    data.append('LONGITUDE_EAST   ' + str(float(e[2])) + '\n')
    data.append('YEAR             ' + str(int(e[4])) + '\n')
    data.append('MONTH            ' + str(int(e[5])) + '\n')
    data.append('DAY              ' + str(int(e[6])) + '\n')
    data.append('HOUR             ' + str(int(e[7])) + '\n')
    data.append('MIN              ' + str(int(e[8])) + '\n')
    data.append('SEC              ' + str(int(e[9])) + '\n')
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
    data.append('BP_F2_STEP1      ' + str(f_c[1]) + '\n')
    data.append('BP_F3_STEP1      ' + str(f_c[2]) + '\n')
    data.append('BP_F4_STEP1      ' + str(f_c[3]) + '\n')
    data.append('DEPTH_UPPERLIM   ' + str(depth_upperlim) + '\n')
    data.append('DEPTH_BOTTOMLIM   ' + str(depth_bottomlim) + '\n')
    data.append('INV_MODE_STEP1    ' + inv_mode_step1 + '\n')
    data.append('LOOPS_SDS_CONF    ' + str(loops_sds_conf) + '\n')
    data.append('PHASES_2_USE_ST1 ' + phases_to_use_st1 + '\n')
    data.append('WEIGHT_A_ST1     ' + str(weight_a_st1) + '\n')
    data.append('WIN_START_A_ST1  ' + str(win_start_a_st1) + '\n')

    with open(inp_file,'w') as file:
        data = file.writelines(data)

    return eve_id,inp_file

def write_kiwi_station_file(station_file,output_file):
    with open(station_file) as f:
        lines = f.readlines()
    data = []
    for i,line in enumerate(lines):
        sta = line.split(' ')
        data.append(str(int(i)) + ' ' +  sta[0] + ' ' + sta[1] + \
                    ' ' + sta[2] + '\n')       
    with open(output_file,'w') as file:
        data = file.writelines(data)    


    
def todisplacement(tr,bitrate,freq4):
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

def parse_kiwi_output(kiwi_output_file,kiwi_bootstrap_file):
    with open(kiwi_output_file) as f:
        lines = f.readlines()
    Mw = float(lines[23].split(' ')[-1])
    mis_fit = float(lines[29].split(' ')[-1])
    no_traces = int(lines[39].split(' ')[6])
    with open(kiwi_bootstrap_file) as f2:
        lines2 = f2.readlines()
    str1 = lines2[2].split(' ')
    M0 = float(str1[2])
    M0_lower = str1[3].split(',')[0]
    M0_lower = float(M0_lower.split('[')[1])
    M0_upper = float(str1[4].split(']')[0])
    M0_error = [M0_lower,M0_upper]
    return Mw,mis_fit,no_traces,M0,M0_error        
    
# main
cat = np.loadtxt(cat_file)
##cat = [cat[8,:]]
id_eqs = []
Mw_ref = []
Ml_ga = []
Ms_ga = []
Mw_sim = []
M0_sim = []
M0_err_sim = []
misfit_sim = []
for e in cat:
    # write the input file for rapidinv12.py
    eve_id,inv_inp_file = write_kiwi_input_file(e,kiwi_work_dir)
##    # down-sample the disp waveforms and save as mseed
##    disp_folder = kiwi_data_dir + eve_id + '/'
##    if not os.path.exists(disp_folder):
##        os.makedirs(disp_folder) 
##    waveform_file = data_dir + 'waveforms_eve_' + str(int(e[0])) + '.pkl'
##    st = read(waveform_file)
##    for tr in st:
##        disp = todisplacement(tr,bitrate,freq4)
##        disp.write(disp_folder+'DISPL.'+ \
##                 tr.stats.station+'.'+tr.stats.channel,format='MSEED')
##    # write the stations.dat file for the event
##    station_file_in = data_dir + 'stations_eve_' + str(int(e[0])) + '.txt'
##    station_file_out = disp_folder + 'stations.dat'
##    write_kiwi_station_file(station_file_in,station_file_out)
##    # run the inversion
##    os.chdir(kiwi_work_dir)
##    os.system('python rapidinv12.py' + ' ' + inv_inp_file)
##    os.chdir(current_dir)
##    pdb.set_trace()
    # parse the kiwi outpit file
    kiwi_output_file = kiwi_result_dir + eve_id + '/step1.earthquakeinfo.dat'
    kiwi_bootstrap_file = kiwi_result_dir + eve_id + '/bootstrap.dat'
    Mw_eve,misfit_eve,no_tr_eve,M0_eve,M0_error_eve = parse_kiwi_output(kiwi_output_file,\
                                                    kiwi_bootstrap_file)
    # generate lists for simulation/reference values
    id_eqs.append(eve_id)
    Mw_ref.append(e[11])
    Ml_ga.append(e[12])
    Ms_ga.append(e[13])
    Mw_sim.append(Mw_eve)
    M0_sim.append(M0_eve)
    M0_err_sim.append(M0_error_eve)
    misfit_sim.append(misfit_eve)
    
# PLOTS
# plot 1: Res: Mw_ref-Mw_sim
res1 = map(sub,Mw_ref,Mw_sim)
res1 = [round(x,2) for x in res1]
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(3,1,1)
xlabels = id_eqs
ax1.plot(range(len(res1)),res1,'o',markersize=10.0,alpha=0.5)
ax1.set_xticks(range(len(res1)))
ax1.set_xticklabels(xlabels,fontsize=12)
ax1.set_xlim([-0.5,len(res1)-0.5])
ax1.plot([-0.5,len(res1)-0.5],[0,0],'k-')
ax1.set_ylim([-0.5,0.5])
ax1.set_yticks(np.arange(-0.5,0.6,0.1))
ax1.grid()
ax1.legend(['This study'])
ax1.set_ylabel('Error : Mw(ref)-Mw(sim)') 

# plot 2: bar plot
x = [0.1,0.2,0.3,0.4]
y1 = len(np.where(np.abs(res1)<=0.1)[0])
y2 = len(np.where((np.abs(res1)<=0.2)&(np.abs(res1)>0.1))[0])
y3 = len(np.where((np.abs(res1)<=0.3)&(np.abs(res1)>0.2))[0])
y4 = len(np.where((np.abs(res1)>0.3))[0])
y_this = np.array([y1,y2,y3,y4])

ax2 = fig.add_subplot(3,1,2)
ax2.bar(x,y_this,facecolor='#9999ff',edgecolor='white',width=0.05,alpha=0.5)
for x1,y1 in zip(x,y_this):
    ax2.text(x1+0.05,y1+0.05,'%d'%y1)

ax2.set_xticks(np.array(x)+0.025)
ax2.set_xticklabels(('<0.1','0.1-0.2','0.2-0.3','>0.3'),fontsize=12)
ax2.set_xlim([0.05,0.5])
ax2.set_ylim([np.min(y_this)-1,np.max(y_this)+1])
ax2.set_ylabel('Frequency')
# plot 4: GOF leveles
ax3 = fig.add_subplot(3,1,3)
xlabels = id_eqs
ax3.plot(range(len(misfit_sim)),misfit_sim,'ko')

patches = []
rect1 = mpatch.Rectangle((-0.5,0),len(xlabels),0.4)
rect2 = mpatch.Rectangle((-0.5,0.4),len(xlabels),0.1)
rect3 = mpatch.Rectangle((-0.5,0.5),len(xlabels),0.1)
rect4 = mpatch.Rectangle((-0.5,0.6),len(xlabels),0.4)
patches.append(rect1)
patches.append(rect2)
patches.append(rect3)
patches.append(rect4)
color = ['green','blue','yellow','red']
collection = PatchCollection(patches,facecolors=color, alpha=0.3)
ax3.add_collection(collection)

ax3.text(len(misfit_sim)-0.6,0.35,'A',fontsize=15, color='red')
ax3.text(len(misfit_sim)-0.6,0.45,'B',fontsize=15, color='red')
ax3.text(len(misfit_sim)-0.6,0.55,'C',fontsize=15, color='red')
ax3.text(len(misfit_sim)-0.6,0.65,'D',fontsize=15, color='red')

ax3.set_xticks(range(len(misfit_sim)))
ax3.set_xticklabels(xlabels,fontsize=12)
ax3.set_xlim([-0.5,len(misfit_sim)-0.5])
ax3.set_ylabel('Goodness-of-fit')
plt.show()


os.system('rm -r ../DATA/auto_*')

