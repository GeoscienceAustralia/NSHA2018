# import modules
import os
import sys
##sys.path.append(os.getcwd())
from modules import *
from obspy.core import UTCDateTime
from obspy.core.event import read_events

# output parameters
# catalogue
cat_out_dir = '../meta-data/catalogue/'
# data
sta_out_dir = '../data/sta_inv/'
wf_out_dir = '../data/wf/'
plot_out_dir = '../plots/'

results_out_dir = '../results/'

# input files
# catalogue
cat_ref_file = '../meta-data/catalogue/cat-ref.csv'
au_shape_file = '../meta-data/aus_shp/AUS_adm0.shp'
# inversion (kiwi)
kiwi_dir = {'kiwi_work_dir':'/home/geoscience/Work/AU_Mw/WORK/',
            'kiwi_data_dir':'/home/geoscience/Work/AU_Mw/DATA/',
            'kiwi_gfdb_dir':'/home/geoscience/Work/AU_Mw/GFDB/',
            'kiwi_result_dir':'/home/geoscience/Work/AU_Mw/RESULTS/'}


# input parameters
# catalogue
search_params_eve_iris = {'start_time':UTCDateTime("2005-01-01T00:00:00"),
                          'end_time':UTCDateTime("2016-01-01T00:00:00"),
                          'min_lat':-44.0, 'max_lat':-9.0,
                          'min_lon':111.0, 'max_lon':156.0,
                          'min_mag':3.5, 'order':'time'}
# data
search_params_data_iris = {'network':"AU,S",'station':"*",'channel':"BH*",
                           't_before_ot':10.*60,'t_after_ot':10.*60, 'maxradius': 10., 
                           'output_dir_sta':sta_out_dir,'output_dir_wf':wf_out_dir}

proc_params = {'freq4':[0.01,0.03,15,20.],'bitrate':2.0,'ratio_vel':20.,'ratio_disp':8.0,
               'freq_range_snr':[0.03,0.1],'per_snr':0.75}


              

# main code
##############################catalogue###########################################

cat_ref = read_cat_ref(cat_ref_file)
##cat_ref.write(cat_out_dir + 'cat_ref.xml', format='QUAKEML')

##cat_iris = get_cat_iris(search_params_eve_iris)
##cat_iris.write(cat_out_dir + 'cat_iris.xml', format='QUAKEML')

##au_poly = get_au_poly(au_shape_file)

##cat_iris_au = filter_cat_poly(cat_iris,au_poly)
##cat_iris_au.write(cat_out_dir + 'cat_iris_au.xml', format='QUAKEML')

##cat1 = read_events(cat_out_dir+'cat_ref.xml')
##cat2 = read_events(cat_out_dir+'cat_iris_au.xml')
##cat_com = merge_cats(cat1,cat2,dt=5.0,dl=10000.)
##cat_com.write(cat_out_dir + 'cat_com_ref_iris.xml', format='QUAKEML')

##############################preparation#########################################
##cat = read_events(cat_out_dir+'cat_com_ref_iris.xml')
cat = cat_ref
cat = cat[7:8]
##get_sta_wf_iris(search_params_data_iris,cat)



##process_wf_cat(cat,sta_out_dir,wf_out_dir,proc_params,plot_out_dir)
output_pre = 'test_'
run_inv_cat(cat,wf_out_dir,sta_out_dir,proc_params,kiwi_dir,output_pre,results_out_dir)


