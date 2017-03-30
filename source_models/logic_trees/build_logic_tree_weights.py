"""Calculate logic tree branch weights based on expert elicitation
workshop - uses expert weights and expert responses

Jonathan Griffin
Geoscience Australia
March 2017
"""

import os, sys
from os.path import join
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

num_ssc_experts = 15
num_gm_experts = 10

calib_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Excalibur/data/'
ssc_weights_file = join(calib_path, 'source.calib.files', 'source.model.calib.weight') #'sourceCalibWeights_format7.m')

ssc_weights = scipy.io.loadmat(ssc_weights_file)
ssc_weights = np.array(ssc_weights['calibration_weight']).flatten()
#print ssc_weights

# Load in target responses
target_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Target_questions'
fig_path =join(target_path, 'seismic_source_results', 'Figures')
ssc_weights
ssc_responses = {}
for i in range(num_ssc_experts):
    print 'Reading response from expert', i+1
    filename = 'Results_expert%i.csv' % (i+1)
    filepath = join(target_path, 'seismic_source_results', filename)
    try:
        data = np.genfromtxt(filepath, delimiter = ',', skip_header = 1, usecols = (1,2,3))#, filling_values = 0, missing_values = 0)
    except IOError:
        print 'Missing file from expert %i!!!' % (i+1)
        continue
    # Check for missing responses
    if np.isnan(data).any():
        print 'Missing values for expert  %i!!!' % (i+1)
        print 'Currently setting all nan to zero'
        data = np.nan_to_num(data)
    weighted_data = data*ssc_weights[i]
    ssc_responses[str(i+1)] = weighted_data
weighted_sum = sum(ssc_responses.values())
#print weighted_sum

# get question ids
ids = []
f_in = open(filepath, 'r')
header = f_in.readline()
for line in f_in.readlines():
    ids.append(line.split(',')[0])
print ids

def get_weights(q_list, weighted_sum):
    ind = []
    for q in q_list:
        index = ids.index(q)
        ind.append(index)
    weights =  weighted_sum[ind][:,1]
    return weights

def autolabel(rects, ax, fontsize = 12):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height+0.05,
                '%.3f' % height,
                ha='center', va='bottom', fontsize = fontsize)

def bar_plot(q_list, weight_list, label_list, filename, title, fontsize = 12, fig_path = fig_path):
    width = 0.35
    x_vals = np.arange(len(q_list))
    fig, ax = plt.subplots()
    rects1 = ax.bar(x_vals, weight_list, width, color='b')
    
    ax.set_ylabel('Weight')
    ax.set_title(title)
    ax.set_xticks(x_vals + width / 2)
    ax.set_xticklabels(label_list)
    ax.set_ylim([0,1])
    ax.set_xlim([-0.5, (len(q_list))])
    autolabel(rects1, ax,fontsize = fontsize )
   
    plt.savefig(join(fig_path, filename))

def bar_subplots(q_lists, weight_lists, label_lists, filename, 
                 title_list, num_row = 2, num_col = 2, 
                 fontsize = 12, fig_path = fig_path):
    """Plot multiple bar plots on one figure
    """
    width = 0.35

    fig = plt.figure()
    for i,weight_list in enumerate(weight_lists):
        x_vals = np.arange(len(q_lists[i]))
        ax = fig.add_subplot(num_row,num_col,i)
        rects1 = ax.bar(x_vals, weight_list, width, color='b')
        ax.set_ylabel('Weight', fontsize = fontsize + 2)
        ax.set_title(title_list[i])
        ax.set_xticks(x_vals + width / 2)
        ax.set_xticklabels(label_lists[i], fontsize = fontsize)
        ax.tick_params(axis='both', which='major', size=fontsize)
        ax.tick_params(axis='both', which='minor', size=fontsize)
        ax.set_ylim([0,1])
        ax.set_xlim([-0.5, (len(q_lists[i]))])
        autolabel(rects1, ax, fontsize = fontsize)
    plt.tight_layout()
    plt.savefig(join(fig_path, filename))

# Source type
src_type =['S2Q1', 'S2Q2', 'S2Q3', 'S2Q4', 'S2Q5']
src_labels = ['Smoothed', 'Background', 'Regional', 'Seismotec', 'Smoothed \n+ faults']
src_type_w =  get_weights(src_type, weighted_sum)
print src_type_w, sum(src_type_w)
bar_plot(src_type, src_type_w, src_labels, 'source_type_weights.png', 'Source Model Type Weights')

# Smoothed seismicity models
ss_models = ['S3Q1', 'S3Q2', 'S3Q3', 'S3Q4']
ss_labels = ['Hall', 'GA fixed', 'GA adaptive', 'Cuthbertson']
ss_w =  get_weights(ss_models, weighted_sum)
print ss_w, sum(ss_w)
bar_plot(ss_models, ss_w, ss_labels, 'smoothed_seismicity_weights.png', 'Smoothed Seismicity Weights')

# Background models
bg_models = ['S4Q1', 'S4Q2', 'S4Q3', 'S4Q4', 'S4Q5']
bg_labels = ['GA NSHA13', 'S&Mc', 'Leonard08', 'Arup', 'Neotectonic \nDomains']
bg_w =  get_weights(bg_models, weighted_sum)
print bg_w, sum(bg_w)
bar_plot(bg_models, bg_w, bg_labels, 'background_model_weights.png', 'Background Model Weights')

# Regional models
rg_models = ['S5Q1', 'S5Q2', 'S5Q3', 'S5Q4']
rg_labels = ['GA NSHA13 \nno H/S', 'GA NSHA13 \nwith H/S', 'AUS6', 'DIM-AUS']
rg_w =  get_weights(rg_models, weighted_sum)
print rg_w, sum(rg_w)
bar_plot(rg_models, rg_w, rg_labels, 'regional_model_weights.png', 'Regional Model Weights')

# Seismotectonic models
st_models = ['S6Q1', 'S6Q2', 'S6Q3']
st_labels = ['GA', 'AUS6', 'DIM-AUS']
st_w =  get_weights(st_models, weighted_sum)
print st_w, sum(st_w)
bar_plot(st_models, st_w, st_labels, 'seismotectonic_model_weights.png', 'Seismotectonic Model Weights')

# Maximum magnitudes
arc_mmaxs = ['S7Q1', 'S7Q2', 'S7Q3', 'S7Q4', 'S7Q5']
arc_mmax_labels = ['7.1', '7.2', '7.3', '7.4', '7.5']
arc_mmax_w = get_weights(arc_mmaxs, weighted_sum)
print arc_mmax_w, sum(arc_mmax_w)
pro_mmaxs = ['S7Q6', 'S7Q7', 'S7Q8', 'S7Q9', 'S7Q10']
pro_mmax_labels = ['7.3', '7.4', '7.5', '7.6', '7.7']
pro_mmax_w = get_weights(pro_mmaxs, weighted_sum)
print pro_mmax_w, sum(pro_mmax_w)
nc_mmaxs = ['S7Q11', 'S7Q12', 'S7Q13', 'S7Q14', 'S7Q15']
nc_mmax_labels = ['7.4', '7.5', '7.6', '7.7', '7.8']
nc_mmax_w = get_weights(nc_mmaxs, weighted_sum)
print nc_mmax_w, sum(nc_mmax_w)
ex_mmaxs = ['S7Q16', 'S7Q17', 'S7Q18', 'S7Q19', 'S7Q20']
ex_mmax_labels = ['7.5', '7.6', '7.7', '7.8', '7.9']
ex_mmax_w = get_weights(ex_mmaxs, weighted_sum)
print ex_mmax_w, sum(ex_mmax_w)

mmax_qlists = [arc_mmaxs, pro_mmaxs, nc_mmaxs, ex_mmaxs]
mmax_weight_lists = [arc_mmax_w, pro_mmax_w, nc_mmax_w, ex_mmax_w]
mmax_label_lists = [arc_mmax_labels, pro_mmax_labels, nc_mmax_labels, ex_mmax_labels]
mmax_title_list = ['Archaen Mmax', 'Proterozoic Mmax', 'Non-cratonic Mmax', 'Extended Mmax']
bar_subplots(mmax_qlists, mmax_weight_lists, mmax_label_lists, 'mmax_weights.png', mmax_title_list)

# b-values
b_qlist = ['S8Q1', 'S8Q2', 'S8Q3', 'S8Q4', 'S8Q5']
b_labels = ['Source zone', 'Neotectonic \ndomain', 'Neotectonic \nsuperdomain', \
            'Continent', 'Model \npropenent']
b_w = get_weights(b_qlist, weighted_sum)
print b_w, sum(b_w)
bar_plot(b_qlist, b_w, b_labels, 'b_value_weights.png', 'b Value Weights')

# Magnitude completeness
mcomp_qlist = ['S9Q1', 'S9Q2', 'S9Q3', 'S9Q4']
mcomp_labels = ['GA NSHA13', 'Cuthbertson', 'Dimas', 'Mote']
mcomp_w = get_weights(mcomp_qlist, weighted_sum)
print mcomp_w, sum(mcomp_w)
bar_plot(mcomp_qlist, mcomp_w, mcomp_labels, \
         'magnitude_completeness_weights.png', \
         'Magnitude Completeness Weights')

# Declustering
ss_dec_qlist = ['S10Q1', 'S10Q2']
ss_dec_labels = ['GA NSHA13', 'Full Catalogue']
ss_dec_w = get_weights(ss_dec_qlist, weighted_sum)
print ss_dec_w, sum(ss_dec_w)
dec_qlist = ['S10Q3', 'S10Q4']
dec_labels = ['GA NSHA13', 'Full Catalogue']
dec_w = get_weights(dec_qlist, weighted_sum)
print dec_w, sum(dec_w)

dec_qlists = [ss_dec_qlist, dec_qlist]
dec_weight_lists = [ss_dec_w, dec_w]
dec_label_lists = [ss_dec_labels, dec_labels]
dec_title_list = ['Smoothed Seismicity Declustering', 'Source Zone Declustering']
bar_subplots(dec_qlists, dec_weight_lists, dec_label_lists, \
             'declustering_weights.png',dec_title_list, \
             num_row=1, num_col=2)

# Fault clustering
cl_c_qlist = ['S1Q1', 'S1Q2']
cl_c_labels = ['Poisson', 'Clustered']
cl_c_w = get_weights(cl_c_qlist, weighted_sum)
print cl_c_w, sum(cl_c_w)
cl_nc_qlist = ['S1Q3', 'S1Q4']
cl_nc_labels = ['Poisson', 'Clustered']
cl_nc_w = get_weights(cl_nc_qlist, weighted_sum)
print cl_nc_w, sum(cl_nc_w)
cl_ex_qlist = ['S1Q5', 'S1Q6']
cl_ex_labels = ['Poisson', 'Clustered']
cl_ex_w = get_weights(cl_ex_qlist, weighted_sum)
print cl_ex_w, sum(cl_ex_w)
cl_type_qlist = ['S1Q7', 'S1Q8']
cl_type_labels = ['Time\n independent', 'Time\n dependent']
cl_type_w = get_weights(cl_type_qlist, weighted_sum)
print cl_type_w, sum(cl_type_w)

cl_qlists = [cl_c_qlist, cl_nc_qlist, cl_ex_qlist, cl_type_qlist]
cl_weight_lists = [cl_c_w, cl_nc_w, cl_ex_w, cl_type_w]
cl_label_lists = [cl_c_labels, cl_nc_labels, cl_ex_labels, cl_type_labels]
cl_title_list = ['Cratonic', 'Non-cratonic', 'Extended', 'Clustering Method']
bar_subplots(cl_qlists, cl_weight_lists, cl_label_lists, 'clustering_weights.png', cl_title_list)

# MFD and integration method
# Cratonic
mfd_c_qlist = ['S1Q9', 'S1Q10', 'S1Q11', 'S1Q12', 'S1Q13', 'S1Q14', 'S1Q15', 'S1Q16', 'S1Q17']
mfd_c_labels = ['YC \nAdd', 'YC \nMB', 'YC \nGeom', 
                'GR \nAdd', 'GR \nMB', 'GR \nGeom',
                'MM \nAdd', 'MM \nMB', 'MM \nGeom']
mfd_c_w = get_weights(mfd_c_qlist, weighted_sum)
print mfd_c_w, sum(mfd_c_w)
# Non-cratonic
mfd_nc_qlist = ['S1Q18', 'S1Q19', 'S1Q20', 'S1Q21', 'S1Q22', 'S1Q23', 'S1Q24', 'S1Q25', 'S1Q26']
mfd_nc_labels = ['YC \nAdd', 'YC \nMB', 'YC \nGeom', 
                'GR \nAdd', 'GR \nMB', 'GR \nGeom',
                'MM \nAdd', 'MM \nMB', 'MM \nGeom']
mfd_nc_w = get_weights(mfd_nc_qlist, weighted_sum)
print mfd_nc_w, sum(mfd_nc_w)
# Extended
mfd_ex_qlist = ['S1Q27', 'S1Q28', 'S1Q29', 'S1Q30', 'S1Q31', 'S1Q32', 'S1Q33', 'S1Q34', 'S1Q35']
mfd_ex_labels = ['YC \nAdd', 'YC \nMB', 'YC \nGeom', 
                'GR \nAdd', 'GR \nMB', 'GR \nGeom',
                'MM \nAdd', 'MM \nMB', 'MM \nGeom']
mfd_ex_w = get_weights(mfd_ex_qlist, weighted_sum)
print mfd_ex_w, sum(mfd_ex_w)

mfd_qlists = [mfd_c_qlist, mfd_nc_qlist, mfd_ex_qlist]
mfd_weight_lists = [mfd_c_w, mfd_nc_w, mfd_ex_w]
mfd_label_lists = [mfd_c_labels, mfd_nc_labels, mfd_ex_labels]
mfd_title_list = ['Cratonic', 'Non-cratonic', 'Extended']
bar_subplots(mfd_qlists, mfd_weight_lists, mfd_label_lists, \
             'fault_mfd_integration_weights.png', mfd_title_list, \
             num_row = 3, num_col = 1, fontsize = 10)
