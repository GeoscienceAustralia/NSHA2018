"""Calculate logic tree branch weights based on expert elicitation
workshop - uses expert weights and expert responses

Jonathan Griffin
Geoscience Australia
March 2017
"""

import os, sys
from os.path import join
import copy
import numpy as np
import scipy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

num_ssc_experts = 15
num_gmm_experts = 10

cw = '0p4'
fig_cw = 'p0.' + cw[-1]

#########################################################
# Seismic source model
calib_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Excalibur/data/'
ssc_weights_file = join(calib_path, 'source.calib.files', 'source.model.calib.weight.cw%s.mat') %cw #'sourceCalibWeights_format7.m')

ssc_weights = scipy.io.loadmat(ssc_weights_file)
#print ssc_weights
ssc_weights = np.array(ssc_weights['source_calibration_weight']).flatten()
#print ssc_weights

# Load in target responses for seismic source mode,l
target_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Target_questions'
fig_path =join(target_path, 'seismic_source_results', 'Figures_%s') % fig_cw
ssc_weights
ssc_responses = {}


for i in range(num_ssc_experts):
    filename = 'Results_expert%i.csv' % (i+1)
    filepath = join(target_path, 'seismic_source_results', filename)
    if i == 0:
        # get question ids
        ids = []
        f_in = open(filepath, 'r')
        header = f_in.readline()
        for line in f_in.readlines():
            ids.append(line.split(',')[0])
        f_in.close()
        S1a = ids[0:2]
        S1b = ids[2:4]
        S1c = ids[4:6]
        S1d = ids[6:8]
        S1e = ids[8:17]
        S1f = ids[17:26]
        S1g = ids[26:35]
        S2 = ids[35:40]
        S3 = ids[40:44]
        S4 = ids[44:49]
        S5 = ids[49:53]
        S6 = ids[53:56]
        S7a = ids[56:61]
        S7b = ids[61:66]
        S7c = ids[66:71]
        S7d = ids[71:76]
        S8 = ids[76:81]
        S9 = ids[81:85]
        S10 = ids[85:89]
        set_list = [S1a,S1b,S1c,S1d,S1e,S1f,S1g,S2,
                    S3,S4,S5,S6,S7a,S7b,S7c,S7d,S8,S9,S10]
        # link ids to index
        set_indices = []
        question_index_list = 0
        for q_set in set_list:
            q_set_indices = []
            for question in q_set:
                q_set_indices.append(question_index_list)
                question_index_list +=1
            set_indices.append(q_set_indices)
        #print set_indices
        # After checking data digitisation, these are known errors where 
        # expert responses don't sum to 1. These will be rescaled
        known_errors = [[1, S4],[4, S4], [6, S7d],[12, S1e], [12, S1f],
                        [12, S1g],[12, S4],[13, S1g],[13, S8],[14,S8]]
    # Read expert responses
    print 'Reading response from expert', i+1
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
   # print data[:,1]
    print 'Sum of weights for expert %i is %.3f' % (i+1, sum(data[:,1]))
    S1_suma = sum(data[0:2,1])
    S1_sumb = sum(data[2:4,1])
    S1_sumc = sum(data[4:6,1])
    S1_sumd = sum(data[6:8,1])
    S1_sume = sum(data[8:17,1])
    S1_sumf = sum(data[17:26,1])
    S1_sumg = sum(data[26:35,1])
    S2_sum = sum(data[35:40,1])
    S3_sum = sum(data[40:44,1])
    S4_sum = sum(data[44:49,1])
    S5_sum = sum(data[49:53,1])
    S6_sum = sum(data[53:56,1])
    S7_suma = sum(data[56:61,1])
    S7_sumb = sum(data[61:66,1])
    S7_sumc = sum(data[66:71,1])
    S7_sumd = sum(data[71:76,1])
    S8_sum = sum(data[76:81,1])
    S9_sum = sum(data[81:85,1])
    S10_sum = sum(data[85:89,1])
    qlist_sum = [S1_suma,S1_sumb,S1_sumc,S1_sumd,S1_sume,S1_sumf,S1_sumg,
                 S2_sum,S3_sum,S4_sum,S5_sum,S6_sum,S7_suma,S7_sumb,S7_sumc,S7_sumd,
                 S8_sum,S9_sum,S10_sum]
    expected_qlist_sum = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 2.]
    for j in range(len(qlist_sum)):
#        print (j+1), qlist_sum[j]
        if np.allclose(qlist_sum[j],expected_qlist_sum[j], rtol=1e-05):
            pass
        else:
#            print (j+1), qlist_sum[j]
            print 'Responses for Question set do not sum to expected %.f' %  expected_qlist_sum[j]
#            print set_list[(j)]
            for ke in known_errors:
                if ke[0] == i+1:
                    if ke[1] == set_list[j]:
                       # print data
                        print 'Known error for expert, normalising values'
                       # print (i+1), set_list[j]
                        original_values = data[set_indices[j][0]:(set_indices[j][-1]+1)]
                        normalised_values = original_values/qlist_sum[j]*expected_qlist_sum[j]
                      #  print original_values, sum(original_values)
                       # print normalised_values, sum(normalised_values)
                        data[set_indices[j][0]:(set_indices[j][-1]+1)] = normalised_values
                       # print data
    weighted_data = data*ssc_weights[i]
    ssc_responses[str(i+1)] = weighted_data
weighted_sum = sum(ssc_responses.values())
#print weighted_sum

######################################################
# Utility functions

def get_weights(q_list, weighted_sum, percentile = None):
    ind = []
    for q in q_list:
        index = ids.index(q)
        ind.append(index)
    weights =  weighted_sum[ind][:,1]
    sum_weights = sum(weights)
    if percentile is not None: 
        # remove low scoring models and re-weight class
        w_sum = 0
        orig_weights = copy.copy(weights)
        new_weight_list = np.zeros(len(orig_weights)) # list to store models we are keeping
        for i in range(len(orig_weights)):
            if w_sum < percentile:
                max_weight_index = np.argmax(orig_weights)
                max_weight = orig_weights[max_weight_index]
                w_sum += max_weight
                new_weight_list[max_weight_index] =  orig_weights[max_weight_index]
                orig_weights[max_weight_index] = 0
            else:
                # Normalise remaining weights
                weights = new_weight_list/sum(new_weight_list)*sum_weights
                break    
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

def bar_plot(q_list, weight_list, label_list, filename, title, fontsize = 12,
             xlabelfont = None, fig_path = fig_path, wide=False, colour_list = None):
    width = 0.35
    x_vals = np.arange(len(q_list))
    fig, ax = plt.subplots()
    if colour_list is not None:
        rects1 = ax.bar(x_vals, weight_list, width, color=colour_list)
    else:
        rects1 = ax.bar(x_vals, weight_list, width, color='b')
    if xlabelfont is None:
        xlabelfont = fontsize
    ax.set_ylabel('Weight')
    ax.set_title(title)
    ax.set_xticks(x_vals + width / 2)
    ax.set_xticklabels(label_list, fontsize=xlabelfont)
    ax.set_ylim([0,1])
    ax.set_xlim([-0.5, (len(q_list))])
    autolabel(rects1, ax,fontsize = fontsize)
    if wide:
        width, height = fig.get_size_inches()
        fig.set_size_inches(width*wide, height)
    plt.tight_layout()
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

def largest_remainder(weights, expected_sum=1, precision=0):
    """Use largest remainder method to round weights such that
    the sum to 1.
    params weights: The raw weights
    params expected_sum: The total the weights should sum to
    params precision: Number of decimal places to round to
    """
    total_number = expected_sum*np.power(10,precision)
    weights = weights*np.power(10,precision)
#    print weights, sum(weights)
    updated_weights = np.floor(weights)
#    print updated_weights, sum(updated_weights)
    remainders = weights - updated_weights
#    print remainders
    unallocated_places = total_number - np.sum(updated_weights)
#    print unallocated_places
    for i in range(int(unallocated_places)):
        max_remainder_index = np.argmax(remainders)
        updated_weights[max_remainder_index] = updated_weights[max_remainder_index] + 1
        remainders[max_remainder_index] = 0
    updated_weights = updated_weights/np.power(10, precision)    
#    print updated_weights
    return updated_weights
    
    

######################################################
# Calculate and plot seismic source model logic tree weights

# Source type
src_type =['S2Q1', 'S2Q2', 'S2Q3', 'S2Q4', 'S2Q5']
src_labels = ['Smoothed', 'Background', 'Regional', 'Seismotec', 'Smoothed \n+ faults']
src_type_w =  get_weights(src_type, weighted_sum)
src_type_w = largest_remainder(src_type_w, expected_sum = 1, precision = 3)
print src_type_w, sum(src_type_w)
bar_plot(src_type, src_type_w, src_labels, 'source_type_weights.png', 'Source Model Type Weights')

# Smoothed seismicity models
ss_models = ['S3Q1', 'S3Q2', 'S3Q3', 'S3Q4']
ss_labels = ['Hall', 'GA fixed', 'GA adaptive', 'Cuthbertson']
ss_w =  get_weights(ss_models, weighted_sum)
ss_w = largest_remainder(ss_w, expected_sum = 1, precision = 3)
print ss_w, sum(ss_w)
bar_plot(ss_models, ss_w, ss_labels, 'smoothed_seismicity_weights.png', 'Smoothed Seismicity Weights')

# Background models
bg_models = ['S4Q1', 'S4Q2', 'S4Q3', 'S4Q4', 'S4Q5']
bg_labels = ['GA NSHA13', 'S&Mc', 'Leonard08', 'Arup', 'Neotectonic \nDomains']
bg_w =  get_weights(bg_models, weighted_sum)
bg_w = largest_remainder(bg_w, expected_sum = 1, precision = 3)
print bg_w, sum(bg_w)
bar_plot(bg_models, bg_w, bg_labels, 'background_model_weights.png', 'Background Model Weights')

# Regional models
rg_models = ['S5Q1', 'S5Q2', 'S5Q3', 'S5Q4']
rg_labels = ['GA NSHA13 \nno H/S', 'GA NSHA13 \nwith H/S', 'AUS6', 'DIM-AUS']
rg_w =  get_weights(rg_models, weighted_sum)
rg_w = largest_remainder(rg_w, expected_sum = 1, precision = 3)
print rg_w, sum(rg_w)
bar_plot(rg_models, rg_w, rg_labels, 'regional_model_weights.png', 'Regional Model Weights')

# Seismotectonic models
st_models = ['S6Q1', 'S6Q2', 'S6Q3']
st_labels = ['GA', 'AUS6', 'DIM-AUS']
st_w =  get_weights(st_models, weighted_sum)
st_w = largest_remainder(st_w, expected_sum = 1, precision = 3)
print st_w, sum(st_w)
bar_plot(st_models, st_w, st_labels, 'seismotectonic_model_weights.png', 'Seismotectonic Model Weights')

# Maximum magnitudes
arc_mmaxs = ['S7Q1', 'S7Q2', 'S7Q3', 'S7Q4', 'S7Q5']
arc_mmax_labels = ['7.1', '7.2', '7.3', '7.4', '7.5']
arc_mmax_w = get_weights(arc_mmaxs, weighted_sum)
arc_mmax_w = largest_remainder(arc_mmax_w, expected_sum = 1, precision = 3)
print arc_mmax_w, sum(arc_mmax_w)
pro_mmaxs = ['S7Q6', 'S7Q7', 'S7Q8', 'S7Q9', 'S7Q10']
pro_mmax_labels = ['7.3', '7.4', '7.5', '7.6', '7.7']
pro_mmax_w = get_weights(pro_mmaxs, weighted_sum)
pro_mmax_w = largest_remainder(pro_mmax_w, expected_sum = 1, precision = 3)
print pro_mmax_w, sum(pro_mmax_w)
nc_mmaxs = ['S7Q11', 'S7Q12', 'S7Q13', 'S7Q14', 'S7Q15']
nc_mmax_labels = ['7.4', '7.5', '7.6', '7.7', '7.8']
nc_mmax_w = get_weights(nc_mmaxs, weighted_sum)
nc_mmax_w = largest_remainder(nc_mmax_w, expected_sum = 1, precision = 3)
print nc_mmax_w, sum(nc_mmax_w)
ex_mmaxs = ['S7Q16', 'S7Q17', 'S7Q18', 'S7Q19', 'S7Q20']
ex_mmax_labels = ['7.5', '7.6', '7.7', '7.8', '7.9']
ex_mmax_w = get_weights(ex_mmaxs, weighted_sum)
ex_mmax_w = largest_remainder(ex_mmax_w, expected_sum = 1, precision = 3)
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
b_w = largest_remainder(b_w, expected_sum = 1, precision = 3)
print b_w, sum(b_w)
bar_plot(b_qlist, b_w, b_labels, 'b_value_weights.png', 'b Value Weights')

# Magnitude completeness
mcomp_qlist = ['S9Q1', 'S9Q2', 'S9Q3', 'S9Q4']
mcomp_labels = ['GA NSHA13', 'Cuthbertson', 'Dimas', 'Mote']
mcomp_w = get_weights(mcomp_qlist, weighted_sum)
mcomp_w = largest_remainder(mcomp_w, expected_sum = 1, precision = 3)
print mcomp_w, sum(mcomp_w)
bar_plot(mcomp_qlist, mcomp_w, mcomp_labels, \
         'magnitude_completeness_weights.png', \
         'Magnitude Completeness Weights')

# Declustering
ss_dec_qlist = ['S10Q1', 'S10Q2']
ss_dec_labels = ['GA NSHA13', 'Full Catalogue']
ss_dec_w = get_weights(ss_dec_qlist, weighted_sum)
ss_dec_w = largest_remainder(ss_dec_w, expected_sum = 1, precision = 3)
print ss_dec_w, sum(ss_dec_w)
dec_qlist = ['S10Q3', 'S10Q4']
dec_labels = ['GA NSHA13', 'Full Catalogue']
dec_w = get_weights(dec_qlist, weighted_sum)
dec_w = largest_remainder(dec_w, expected_sum = 1, precision = 3)
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
cl_c_w = largest_remainder(cl_c_w, expected_sum = 1, precision = 3)
print cl_c_w, sum(cl_c_w)
cl_nc_qlist = ['S1Q3', 'S1Q4']
cl_nc_labels = ['Poisson', 'Clustered']
cl_nc_w = get_weights(cl_nc_qlist, weighted_sum)
cl_nc_w = largest_remainder(cl_nc_w, expected_sum = 1, precision = 3)
print cl_nc_w, sum(cl_nc_w)
cl_ex_qlist = ['S1Q5', 'S1Q6']
cl_ex_labels = ['Poisson', 'Clustered']
cl_ex_w = get_weights(cl_ex_qlist, weighted_sum)
cl_ex_w = largest_remainder(cl_ex_w, expected_sum = 1, precision = 3)
print cl_ex_w, sum(cl_ex_w)
cl_type_qlist = ['S1Q7', 'S1Q8']
cl_type_labels = ['Time\n independent', 'Time\n dependent']
cl_type_w = get_weights(cl_type_qlist, weighted_sum)
cl_type_w = largest_remainder(cl_type_w, expected_sum = 1, precision = 3)
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
mfd_c_w = largest_remainder(mfd_c_w, expected_sum = 1, precision = 3)
print mfd_c_w, sum(mfd_c_w)
# Non-cratonic
mfd_nc_qlist = ['S1Q18', 'S1Q19', 'S1Q20', 'S1Q21', 'S1Q22', 'S1Q23', 'S1Q24', 'S1Q25', 'S1Q26']
mfd_nc_labels = ['YC \nAdd', 'YC \nMB', 'YC \nGeom', 
                'GR \nAdd', 'GR \nMB', 'GR \nGeom',
                'MM \nAdd', 'MM \nMB', 'MM \nGeom']
mfd_nc_w = get_weights(mfd_nc_qlist, weighted_sum)
mfd_nc_w = largest_remainder(mfd_nc_w, expected_sum = 1, precision = 3)
print mfd_nc_w, sum(mfd_nc_w)
# Extended
mfd_ex_qlist = ['S1Q27', 'S1Q28', 'S1Q29', 'S1Q30', 'S1Q31', 'S1Q32', 'S1Q33', 'S1Q34', 'S1Q35']
mfd_ex_labels = ['YC \nAdd', 'YC \nMB', 'YC \nGeom', 
                'GR \nAdd', 'GR \nMB', 'GR \nGeom',
                'MM \nAdd', 'MM \nMB', 'MM \nGeom']
mfd_ex_w = get_weights(mfd_ex_qlist, weighted_sum)
mfd_ex_w = largest_remainder(mfd_ex_w, expected_sum = 1, precision = 3)
print mfd_ex_w, sum(mfd_ex_w)

mfd_qlists = [mfd_c_qlist, mfd_nc_qlist, mfd_ex_qlist]
mfd_weight_lists = [mfd_c_w, mfd_nc_w, mfd_ex_w]
mfd_label_lists = [mfd_c_labels, mfd_nc_labels, mfd_ex_labels]
mfd_title_list = ['Cratonic', 'Non-cratonic', 'Extended']
bar_subplots(mfd_qlists, mfd_weight_lists, mfd_label_lists, \
             'fault_mfd_integration_weights.png', mfd_title_list, \
             num_row = 3, num_col = 1, fontsize = 10)

######################################################
# Calculate and plot ground motion logic tree weights
# Load data
calib_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Excalibur/data/GMPE/'
gmm_weights_file = join(calib_path, 'ground.motion.model.calib.weight.%s.mat') % fig_cw #'sourceCalibWeights_format7.m')

gmm_weights = scipy.io.loadmat(gmm_weights_file)
print gmm_weights
gmm_weights = np.array(gmm_weights['gmm_calibration_weight']).flatten()
#print gmm_weights

# Load in target responses for ground motion model
print 'Loading ground motion model weight data'
target_path = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/Expert_Elicitation/Target_questions'
#fig_path =join(target_path, 'ground_motion_results', 'Figures_p0.3')
#gmm_weights
gmm_responses = {}
for i in range(num_gmm_experts):
    filename = 'Results_expert%i.csv' % (i+1)
    filepath = join(target_path, 'ground_motion_results', filename)
    if i == 0:
        # get question ids
        ids = []
        f_in = open(filepath, 'r')
        header = f_in.readline()
        for line in f_in.readlines():
            ids.append(line.split(',')[0])
        f_in.close()
        S1a = ids[0:4]
        S1b = ids[4:8]
        S1c = ids[8:13]
        S2a = ids[13:17]
        S2b = ids[17:24]
        S2c = ids[24:27]
        S2d = ids[27:30]
        S3a = ids[30:34]
        S3b = ids[34:41]
        S3c = ids[41:44]
        S3d = ids[44:47]
        S4a = ids[47:51]
        S4b = ids[51:58]
        S4c = ids[58:61]
        S4d = ids[61:64]
        S4e = ids[64:68]
        S5a = ids[68:73]
        S5b = ids[73:78]
        set_list = [S1a,S1b,S1c,S2a,S2b,S2c,S2d,S3a,S3b,S3c,S3d,
                    S4a,S4b,S4c,S4d,S4e,S5a,S5b]
        for q_set in set_list:
            print q_set
        # link ids to index
        set_indices = []
        question_index_list = 0
        for q_set in set_list:
            q_set_indices = []
            for question in q_set:
                q_set_indices.append(question_index_list)
                question_index_list +=1
            set_indices.append(q_set_indices)
        #print set_indices
        # After checking data digitisation, these are known errors where 
        # expert responses don't sum to 1. These will be rescaled
        known_errors = [[4,S2a], [2, S1b], [7,S4e]]
    print 'Reading response from expert', i+1
    filename = 'Results_expert%i.csv' % (i+1)
    filepath = join(target_path, 'ground_motion_results', filename)
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

    print 'Sum of weights for expert %i is %.3f' % (i+1, sum(data[:,1]))
    S1_suma = sum(data[0:4,1])
    S1_sumb = sum(data[4:8,1])
    S1_sumc = sum(data[8:13,1])
    S2_suma = sum(data[13:17,1])
    S2_sumb = sum(data[17:24,1])
    S2_sumc = sum(data[24:27,1])
    S2_sumd = sum(data[27:30,1])
    S3_suma = sum(data[30:34,1])
    S3_sumb = sum(data[34:41,1])
    S3_sumc = sum(data[41:44,1])
    S3_sumd = sum(data[44:47,1])
    S4_suma = sum(data[47:51,1])
    S4_sumb = sum(data[51:58,1])
    S4_sumc = sum(data[58:61,1])
    S4_sumd = sum(data[61:64,1])
    S4_sume = sum(data[64:68,1])
    S5_suma = sum(data[68:73,1])
    S5_sumb = sum(data[73:78,1])
    qlist_sum = [S1_suma,S1_sumb,S1_sumc,S2_suma,S2_sumb,S2_sumc,S2_sumd,
                 S3_suma,S3_sumb,S3_sumc,S3_sumd,S4_suma,S4_sumb,S4_sumc,
                 S4_sumd,S4_sume,S5_suma,S5_sumb]
    expected_qlist_sum = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
    for j in range(len(qlist_sum)):
#        print (j+1), qlist_sum[j]
        if np.allclose(qlist_sum[j],expected_qlist_sum[j], rtol=1e-05):
            pass
        else:
            print (j+1), qlist_sum[j]
            print 'Responses for Question set do not sum to expected %.f' %  expected_qlist_sum[j]
            print set_list[(j)]
            for ke in known_errors:
                if ke[0] == i+1:
                    if ke[1] == set_list[j]:
                       # print data
                        print 'Known error for expert, normalising values'
                       # print (i+1), set_list[j]
                        original_values = data[set_indices[j][0]:(set_indices[j][-1]+1)]
                        normalised_values = original_values/qlist_sum[j]*expected_qlist_sum[j]
                      #  print original_values, sum(original_values)
                       # print normalised_values, sum(normalised_values)
                        data[set_indices[j][0]:(set_indices[j][-1]+1)] = normalised_values
                       # print data
    weighted_data = data*gmm_weights[i]
    gmm_responses[str(i+1)] = weighted_data
weighted_sum = sum(gmm_responses.values())
#print weighted_sum

# get question ids
ids = []
f_in = open(filepath, 'r')
header = f_in.readline()
for line in f_in.readlines():
    ids.append(line.split(',')[0])
print ids

c_gmm_model_all_w = []
c_gmm_model_all_labels = []
nc_ex_gmm_model_all_w = []
nc_ex_gmm_model_all_labels = []
banda_gmm_model_all_w = []
banda_gmm_model_all_labels = []

gmm_fig_path = fig_path =join(target_path, 'ground_motion_results', 'Figures_%s') % fig_cw
# Cratonic GMM type
c_gmm_region_qlist = ['S1Q1', 'S1Q2', 'S1Q3', 'S1Q4']
c_gmm_region_labels = ['Australian', 'CEUS', 'California', 'European']
c_gmm_region_w = get_weights(c_gmm_region_qlist, weighted_sum)
c_gmm_region_w = largest_remainder(c_gmm_region_w, expected_sum = 1, precision = 3)
print c_gmm_region_w, sum(c_gmm_region_w)
bar_plot(c_gmm_region_qlist, c_gmm_region_w, c_gmm_region_labels, \
         'cratonic_gmm_region_weights.png', \
         'Cratonic GMM Region Weights', fig_path = gmm_fig_path)

# Non-cratonic and extended GMM type
nc_ex_gmm_region_qlist = ['S1Q5', 'S1Q6', 'S1Q7', 'S1Q8']
nc_ex_gmm_region_labels = ['Australian', 'CEUS', 'California', 'European']
nc_ex_gmm_region_w = get_weights(nc_ex_gmm_region_qlist, weighted_sum)
nc_ex_gmm_region_w = largest_remainder(nc_ex_gmm_region_w, expected_sum = 1, precision = 3)
print nc_ex_gmm_region_w, sum(nc_ex_gmm_region_w)
bar_plot(nc_ex_gmm_region_qlist, nc_ex_gmm_region_w, nc_ex_gmm_region_labels, \
         'non_cratonic_extended_gmm_region_weights.png', \
         'Non-cratonic and Extended GMM Region Weights', fig_path = gmm_fig_path)

# Banda Sea GMM type
banda_gmm_region_qlist = ['S1Q9', 'S1Q10', 'S1Q11', 'S1Q12', 'S1Q13']
banda_gmm_region_labels = ['Australian', 'CEUS', 'California', 'European','Subduction \nIntraslab']
banda_gmm_region_w = get_weights(banda_gmm_region_qlist, weighted_sum)
banda_gmm_region_w = largest_remainder(banda_gmm_region_w, expected_sum = 1, precision = 3)
print banda_gmm_region_w, sum(banda_gmm_region_w)
bar_plot(banda_gmm_region_qlist, banda_gmm_region_w, banda_gmm_region_labels, \
         'banda_sea_gmm_region_weights.png', \
         'Banda Sea GMM Region Weights', fig_path = gmm_fig_path)

##############
# Cratonic models
# Cratonic GMM Australian models
c_gmm_aust_qlist = ['S2Q1', 'S2Q2', 'S2Q3', 'S2Q4']
c_gmm_aust_labels = ['Allen2012', 'Allen2012 \n RedSigma', 'Som2009 \n Non-cratonic', 'Som2009 \n Yilgarn']
c_gmm_aust_w = get_weights(c_gmm_aust_qlist, weighted_sum)
# Remove Allen2012 reduced sigma model for now, renormalise other values
c_gmm_aust_w = np.delete(c_gmm_aust_w,1)
del c_gmm_aust_qlist[1]
del c_gmm_aust_labels[1]
c_gmm_aust_w = c_gmm_aust_w/sum(c_gmm_aust_w)
c_gmm_aust_w = largest_remainder(c_gmm_aust_w, expected_sum = 1, precision = 3)
print c_gmm_aust_w, sum(c_gmm_aust_w)
c_gmm_aust_region_w = c_gmm_aust_w *c_gmm_region_w[0]
c_gmm_model_all_w += list(c_gmm_aust_region_w)
c_gmm_model_all_labels += c_gmm_aust_labels
bar_plot(c_gmm_aust_qlist, c_gmm_aust_w, c_gmm_aust_labels, \
         'cratonic_gmm_australian_weights.png', \
         'Cratonic Australian GMM Weights', fig_path = gmm_fig_path)

# Cratonic GMM CEUS models
c_gmm_ceus_qlist = ['S2Q5', 'S2Q6', 'S2Q7', 'S2Q8', 'S2Q9', 'S2Q10', 'S2Q11']
c_gmm_ceus_labels = ['AtkBoore\n2006', 'AtkBoore\n2006\nMod2011', 'Campbell\n2003', 'Pezeshk\n2011', 'Silva2002\nMwNSHMP\n2008', 'Toro\n2002', 'YenAtk\n2015']
c_gmm_ceus_w = get_weights(c_gmm_ceus_qlist, weighted_sum)
c_gmm_ceus_w = largest_remainder(c_gmm_ceus_w, expected_sum = 1, precision = 3)
print c_gmm_ceus_w, sum(c_gmm_ceus_w)
c_gmm_ceus_region_w = c_gmm_ceus_w *c_gmm_region_w[1]
c_gmm_model_all_w += list(c_gmm_ceus_region_w)
c_gmm_model_all_labels += c_gmm_ceus_labels
bar_plot(c_gmm_ceus_qlist, c_gmm_ceus_w, c_gmm_ceus_labels, \
         'cratonic_gmm_CEUS_weights.png', \
         'Cratonic CEUS GMM Weights', fig_path = gmm_fig_path)

# Cratonic GMM California models
c_gmm_cal_qlist = ['S2Q12', 'S2Q13', 'S2Q14']
c_gmm_cal_labels = ['Boore2014', 'ChiYou2008', 'ChiYou2014']
c_gmm_cal_w = get_weights(c_gmm_cal_qlist, weighted_sum)
c_gmm_cal_w = largest_remainder(c_gmm_cal_w, expected_sum = 1, precision = 3)
print c_gmm_cal_w, sum(c_gmm_cal_w)
c_gmm_cal_region_w = c_gmm_cal_w *c_gmm_region_w[2]
c_gmm_model_all_w += list(c_gmm_cal_region_w)
c_gmm_model_all_labels += c_gmm_cal_labels
bar_plot(c_gmm_cal_qlist, c_gmm_cal_w, c_gmm_cal_labels, \
         'cratonic_gmm_california_weights.png', \
         'Cratonic Californian GMM Weights', fig_path = gmm_fig_path)

# Cratonic GMM Eurpoean models
c_gmm_eur_qlist = ['S2Q15', 'S2Q16', 'S2Q17']
c_gmm_eur_labels = ['ChiYo2008\nSWISS01', 'Rietbrock\n2013SS', 'Zhao2006\nAscSWISS05']
c_gmm_eur_w = get_weights(c_gmm_eur_qlist, weighted_sum)
c_gmm_eur_w = largest_remainder(c_gmm_eur_w, expected_sum = 1, precision = 3)
print c_gmm_eur_w, sum(c_gmm_eur_w)
c_gmm_eur_region_w = c_gmm_eur_w *c_gmm_region_w[3]
c_gmm_model_all_w += list(c_gmm_eur_region_w)
c_gmm_model_all_labels += c_gmm_eur_labels
bar_plot(c_gmm_eur_qlist, c_gmm_eur_w, c_gmm_eur_labels, \
         'cratonic_gmm_european_weights.png', \
         'Cratonic Eurpean GMM Weights', fig_path = gmm_fig_path)
# Collate for percentile calcs later
c_gmm_w_list = [c_gmm_aust_region_w, c_gmm_ceus_region_w, c_gmm_cal_region_w, c_gmm_eur_region_w]

##############
# Non-cratonic and Extended models
# Non-cratonic and Extended GMM Australian models
nc_ex_gmm_aust_qlist = ['S3Q1', 'S3Q2', 'S3Q3', 'S3Q4']
nc_ex_gmm_aust_labels = ['Allen2012', 'Allen2012 \n RedSigma', 'Som2009 \n Non-cratonic', 'Som2009 \n Yilgarn']
nc_ex_gmm_aust_w = get_weights(nc_ex_gmm_aust_qlist, weighted_sum)
# Remove Allen2012 reduced sigma model for now, renormalise other values
nc_ex_gmm_aust_w = np.delete(nc_ex_gmm_aust_w, 1)
del nc_ex_gmm_aust_qlist[1]
del nc_ex_gmm_aust_labels[1]
nc_ex_gmm_aust_w = nc_ex_gmm_aust_w/sum(nc_ex_gmm_aust_w)
nc_ex_gmm_aust_w = largest_remainder(nc_ex_gmm_aust_w, expected_sum = 1, precision = 3)
print nc_ex_gmm_aust_w, sum(nc_ex_gmm_aust_w)
nc_ex_gmm_aust_region_w = nc_ex_gmm_aust_w *nc_ex_gmm_region_w[0]
nc_ex_gmm_model_all_w += list(nc_ex_gmm_aust_region_w)
nc_ex_gmm_model_all_labels += nc_ex_gmm_aust_labels
bar_plot(nc_ex_gmm_aust_qlist, nc_ex_gmm_aust_w, nc_ex_gmm_aust_labels, \
         'nc_ex_gmm_australian_weights.png', \
         'Non-cratonic and Extended Australian GMM Weights', fig_path = gmm_fig_path)

# Non-cratonic and Extended GMM CEUS models
nc_ex_gmm_ceus_qlist = ['S3Q5', 'S3Q6', 'S3Q7', 'S3Q8', 'S3Q9', 'S3Q10', 'S3Q11']
nc_ex_gmm_ceus_labels = ['AtkBoore\n2006', 'AtkBoore\n2006\nMod2011', 'Campbell\n2003', 'Pezeshk\n2011', 'Silva2002\nMwNSHMP\n2008', 'Toro\n2002', 'YenAtk\n2015']
nc_ex_gmm_ceus_w = get_weights(nc_ex_gmm_ceus_qlist, weighted_sum)
nc_ex_gmm_ceus_w = largest_remainder(nc_ex_gmm_ceus_w, expected_sum = 1, precision = 3)
print nc_ex_gmm_ceus_w, sum(nc_ex_gmm_ceus_w)
nc_ex_gmm_ceus_region_w = nc_ex_gmm_ceus_w *nc_ex_gmm_region_w[1]
nc_ex_gmm_model_all_w += list(nc_ex_gmm_ceus_region_w)
nc_ex_gmm_model_all_labels += nc_ex_gmm_ceus_labels
bar_plot(nc_ex_gmm_ceus_qlist, nc_ex_gmm_ceus_w, nc_ex_gmm_ceus_labels, \
         'nc_ex_gmm_CEUS_weights.png', \
         'Non-cratonic and Extended CEUS GMM Weights', fig_path = gmm_fig_path)

# Non-cratonic and Extended GMM California models
nc_ex_gmm_cal_qlist = ['S3Q12', 'S3Q13', 'S3Q14']
nc_ex_gmm_cal_labels = ['Boore2014', 'ChiYou2008', 'ChiYou2014']
nc_ex_gmm_cal_w = get_weights(nc_ex_gmm_cal_qlist, weighted_sum)
nc_ex_gmm_cal_w = largest_remainder(nc_ex_gmm_cal_w, expected_sum = 1, precision = 3)
print nc_ex_gmm_cal_w, sum(nc_ex_gmm_cal_w)
nc_ex_gmm_cal_region_w = nc_ex_gmm_cal_w *nc_ex_gmm_region_w[2]
nc_ex_gmm_model_all_w += list(nc_ex_gmm_cal_region_w)
nc_ex_gmm_model_all_labels += nc_ex_gmm_cal_labels
bar_plot(nc_ex_gmm_cal_qlist, nc_ex_gmm_cal_w, nc_ex_gmm_cal_labels, \
         'nc_ex_gmm_california_weights.png', \
         'Non-cratonic and Extended Californian GMM Weights', fig_path = gmm_fig_path)

# Non-cratonic and Extended GMM Eurpoean models
nc_ex_gmm_eur_qlist = ['S3Q15', 'S3Q16', 'S3Q17']
nc_ex_gmm_eur_labels = ['ChiYo2008\nSWISS01', 'Rietbrock\n2013SS', 'Zhao2006\nAscSWISS05']
nc_ex_gmm_eur_w = get_weights(nc_ex_gmm_eur_qlist, weighted_sum)
nc_ex_gmm_eur_w = largest_remainder(nc_ex_gmm_eur_w, expected_sum = 1, precision = 3)
print nc_ex_gmm_eur_w, sum(nc_ex_gmm_eur_w)
nc_ex_gmm_eur_region_w = nc_ex_gmm_eur_w *nc_ex_gmm_region_w[3]
nc_ex_gmm_model_all_w += list(nc_ex_gmm_eur_region_w)
nc_ex_gmm_model_all_labels += nc_ex_gmm_eur_labels
bar_plot(nc_ex_gmm_eur_qlist, nc_ex_gmm_eur_w, nc_ex_gmm_eur_labels, \
         'nc_ex_gmm_european_weights.png', \
         'Non-cratonic and Extended Eurpean GMM Weights', fig_path = gmm_fig_path)
# Collate for percentile calcs later
nc_ex_gmm_w_list = [nc_ex_gmm_aust_region_w, nc_ex_gmm_ceus_region_w, nc_ex_gmm_cal_region_w, nc_ex_gmm_eur_region_w]

##############
# Banda Sea models
# Banda Sea GMM Australian models
banda_gmm_aust_qlist = ['S4Q1', 'S4Q2', 'S4Q3', 'S4Q4']
banda_gmm_aust_labels = ['Allen2012', 'Allen2012 \n RedSigma', 'Som2009 \n Non-cratonic', 'Som2009 \n Yilgarn']
banda_gmm_aust_w = get_weights(banda_gmm_aust_qlist, weighted_sum)
# Remove Allen2012 reduced sigma model for now, renormalise other values
banda_gmm_aust_w = np.delete(banda_gmm_aust_w, 1)
del banda_gmm_aust_qlist[1]
del banda_gmm_aust_labels[1]
banda_gmm_aust_w = banda_gmm_aust_w/sum(banda_gmm_aust_w)
banda_gmm_aust_w = largest_remainder(banda_gmm_aust_w, expected_sum = 1, precision = 3)
print banda_gmm_aust_w, sum(banda_gmm_aust_w)
banda_gmm_aust_region_w = banda_gmm_aust_w *banda_gmm_region_w[0]
banda_gmm_model_all_w += list(banda_gmm_aust_region_w)
banda_gmm_model_all_labels += banda_gmm_aust_labels
bar_plot(banda_gmm_aust_qlist, banda_gmm_aust_w, banda_gmm_aust_labels, \
         'banda_gmm_australian_weights.png', \
         'Banda Sea Australian GMM Weights', fig_path = gmm_fig_path)

# Banda Sea GMM CEUS models
banda_gmm_ceus_qlist = ['S4Q5', 'S4Q6', 'S4Q7', 'S4Q8', 'S4Q9', 'S4Q10', 'S4Q11']
banda_gmm_ceus_labels = ['AtkBoore\n2006', 'AtkBoore\n2006\nMod2011', 'Campbell\n2003', 'Pezeshk\n2011', 'Silva2002\nMwNSHMP\n2008', 'Toro\n2002', 'YenAtk\n2015']
banda_gmm_ceus_w = get_weights(banda_gmm_ceus_qlist, weighted_sum)
banda_gmm_ceus_w = largest_remainder(banda_gmm_ceus_w, expected_sum = 1, precision = 3)
print banda_gmm_ceus_w, sum(banda_gmm_ceus_w)
banda_gmm_ceus_region_w = banda_gmm_ceus_w *banda_gmm_region_w[1]
banda_gmm_model_all_w += list(banda_gmm_ceus_region_w)
banda_gmm_model_all_labels += banda_gmm_ceus_labels
bar_plot(banda_gmm_ceus_qlist, banda_gmm_ceus_w, banda_gmm_ceus_labels, \
         'banda_gmm_CEUS_weights.png', \
         'Banda Sea CEUS GMM Weights', fig_path = gmm_fig_path)

# Banda Sea GMM California models
banda_gmm_cal_qlist = ['S4Q12', 'S4Q13', 'S4Q14']
banda_gmm_cal_labels = ['Boore2014', 'ChiYou2008', 'ChiYou2014']
banda_gmm_cal_w = get_weights(banda_gmm_cal_qlist, weighted_sum)
banda_gmm_cal_w = largest_remainder(banda_gmm_cal_w, expected_sum = 1, precision = 3)
print banda_gmm_cal_w, sum(banda_gmm_cal_w)
banda_gmm_cal_region_w = banda_gmm_cal_w *banda_gmm_region_w[2]
banda_gmm_model_all_w += list(banda_gmm_cal_region_w)
banda_gmm_model_all_labels += banda_gmm_cal_labels
bar_plot(banda_gmm_cal_qlist, banda_gmm_cal_w, banda_gmm_cal_labels, \
         'banda_gmm_california_weights.png', \
         'Banda Sea Californian GMM Weights', fig_path = gmm_fig_path)

# Banda Sea GMM Eurpoean models
banda_gmm_eur_qlist = ['S4Q15', 'S4Q16', 'S4Q17']
banda_gmm_eur_labels = ['ChiYo2008\nSWISS01', 'Rietbrock\n2013SS', 'Zhao2006\nAscSWISS05']
banda_gmm_eur_w = get_weights(banda_gmm_eur_qlist, weighted_sum)
banda_gmm_eur_w = largest_remainder(banda_gmm_eur_w, expected_sum = 1, precision = 3)
print banda_gmm_eur_w, sum(banda_gmm_eur_w)
banda_gmm_eur_region_w = banda_gmm_eur_w *banda_gmm_region_w[3]
banda_gmm_model_all_w += list(banda_gmm_eur_region_w)
banda_gmm_model_all_labels += banda_gmm_eur_labels
bar_plot(banda_gmm_eur_qlist, banda_gmm_eur_w, banda_gmm_eur_labels, \
         'banda_gmm_european_weights.png', \
         'Banda Sea Eurpean GMM Weights', fig_path = gmm_fig_path)

# Banda Sea GMM intraslab models
#percentile = 0.75
banda_gmm_inslab_qlist = ['S4Q18', 'S4Q19', 'S4Q20', 'S4Q21']
banda_gmm_inslab_labels = ['Abrahamson\n2015SSlab', 'AtkBoore\n2003SSlab', 'Garcia\n2005SSlab', 'MegaPan\n2010']
banda_gmm_inslab_w = get_weights(banda_gmm_inslab_qlist, weighted_sum)#, percentile = percentile)
banda_gmm_inslab_w = largest_remainder(banda_gmm_inslab_w, expected_sum = 1, precision = 3)
print banda_gmm_inslab_w, sum(banda_gmm_inslab_w)
banda_gmm_inslab_region_w = banda_gmm_inslab_w *banda_gmm_region_w[4]
banda_gmm_model_all_w += list(banda_gmm_inslab_region_w)
banda_gmm_model_all_labels += banda_gmm_inslab_labels
bar_plot(banda_gmm_inslab_qlist, banda_gmm_inslab_w, banda_gmm_inslab_labels, \
         'banda_gmm_inslab_weights.png', \
         'Banda Sea Intraslab GMM Weights', fig_path = gmm_fig_path)
# Collate for percentile calcs later
banda_gmm_w_list = [banda_gmm_aust_region_w, banda_gmm_ceus_region_w, banda_gmm_cal_region_w, banda_gmm_eur_region_w, banda_gmm_inslab_region_w]

##########################
# GMM Cutoff distances
# Cratonic, Non-Cratonic and Extended cutoff distances
c_nc_ex_cutoff_qlist = ['S5Q1', 'S5Q2', 'S5Q3', 'S5Q4', 'S5Q5']
c_nc_ex_cutoff_labels = ['200', '300', '400', '500', '600']
c_nc_ex_cutoff_w = get_weights(c_nc_ex_cutoff_qlist, weighted_sum)
c_nc_ex_cutoff_w = largest_remainder(c_nc_ex_cutoff_w, expected_sum = 1, precision = 3)
print c_nc_ex_cutoff_w, sum(c_nc_ex_cutoff_w)
bar_plot(c_nc_ex_cutoff_qlist, c_nc_ex_cutoff_w, c_nc_ex_cutoff_labels, \
         'c_nc_ex_cutoff_distance_weights.png', \
         'Cratonic, Non-cratonic and Extended Cutoff Distance Weights', fig_path = gmm_fig_path)

# Banda Sea cutoff distances
banda_cutoff_qlist = ['S5Q6', 'S5Q7', 'S5Q8', 'S5Q9']
banda_cutoff_labels = ['500', '600', '800', '1000']
banda_cutoff_w = get_weights(banda_cutoff_qlist, weighted_sum)
banda_cutoff_w = largest_remainder(banda_cutoff_w, expected_sum = 1, precision = 3)
print banda_cutoff_w, sum(banda_cutoff_w)
bar_plot(banda_cutoff_qlist, banda_cutoff_w, banda_cutoff_labels, \
         'banda_cutoff_distance_weights.png', \
         'Banda Sea Cutoff Distance Weights', fig_path = gmm_fig_path)

#################################
# Plot full model weights 
bar_plot(banda_gmm_model_all_w, banda_gmm_model_all_w, banda_gmm_model_all_labels, \
         'banda_gmm_all_weights.png', \
         'Banda Sea All Model Weights', fig_path = gmm_fig_path, fontsize=10, wide=2.5)
bar_plot(c_gmm_model_all_w, c_gmm_model_all_w, c_gmm_model_all_labels, \
         'cratonic_gmm_all_weights.png', \
         'Cratonic All Model Weights', fig_path = gmm_fig_path, fontsize=10, wide=2)
bar_plot(nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_labels, \
         'noncratonic_extended_gmm_all_weights.png', \
         'Non-cratonic and Extended All Model Weights', fig_path = gmm_fig_path, fontsize=10, wide=2)
################################
# Sort by weight and identy nth percentile cutoff
import copy
percentiles = [0.75, 0.8]
colour_in = 'b'
colour_out = 'r'
for percentile in percentiles:
    print 'Cut-off percentile of', (100*percentile)
    banda_gmm_w_sum = 0
    banda_gmm_model_all_w_copy = copy.copy(banda_gmm_model_all_w)
    banda_gmm_colour_list = [colour_out]*len(banda_gmm_model_all_w)
    banda_gmm_reduced_list = [] # list to store models we are keeping
    banda_gmm_reduced_labels = []
    for i in range(len(banda_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(banda_gmm_model_all_w_copy)
        max_weight = banda_gmm_model_all_w_copy[max_weight_index]
        banda_gmm_w_sum += max_weight
        if banda_gmm_w_sum < percentile:
            # Weights we are keeping
            banda_gmm_colour_list[max_weight_index] = colour_in
            banda_gmm_reduced_list.append(banda_gmm_model_all_w_copy[max_weight_index])
            banda_gmm_reduced_labels.append(banda_gmm_model_all_labels[max_weight_index])
            banda_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(banda_gmm_model_all_w, banda_gmm_model_all_w, banda_gmm_model_all_labels, \
             'banda_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'Banda Sea All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = banda_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    banda_gmm_reduced_list = banda_gmm_reduced_list/sum(banda_gmm_reduced_list)
#    print banda_gmm_reduced_list, sum(banda_gmm_reduced_list)
    bar_plot(banda_gmm_reduced_list, banda_gmm_reduced_list, banda_gmm_reduced_labels, \
             'banda_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'Banda Sea Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)
    #############
    c_gmm_w_sum = 0
    c_gmm_model_all_w_copy = copy.copy(c_gmm_model_all_w)
    c_gmm_colour_list = [colour_out]*len(c_gmm_model_all_w)
    c_gmm_reduced_list = [] # list to store models we are keeping
    c_gmm_reduced_labels = []
    for i in range(len(c_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(c_gmm_model_all_w_copy)
        max_weight = c_gmm_model_all_w_copy[max_weight_index]
        c_gmm_w_sum += max_weight
        if c_gmm_w_sum < percentile:
            # Weights we are keeping
            c_gmm_colour_list[max_weight_index] = colour_in
            c_gmm_reduced_list.append(c_gmm_model_all_w_copy[max_weight_index])
            c_gmm_reduced_labels.append(c_gmm_model_all_labels[max_weight_index])
            c_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(c_gmm_model_all_w, c_gmm_model_all_w, c_gmm_model_all_labels, \
             'c_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'Cratonic All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = c_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    c_gmm_reduced_list = c_gmm_reduced_list/sum(c_gmm_reduced_list)
    bar_plot(c_gmm_reduced_list, c_gmm_reduced_list, c_gmm_reduced_labels, \
             'c_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'Cratonic Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)
    #############
    nc_ex_gmm_w_sum = 0
    nc_ex_gmm_model_all_w_copy = copy.copy(nc_ex_gmm_model_all_w)
    nc_ex_gmm_colour_list = [colour_out]*len(nc_ex_gmm_model_all_w)
    nc_ex_gmm_reduced_list = [] # list to store models we are keeping
    nc_ex_gmm_reduced_labels = []
    for i in range(len(nc_ex_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(nc_ex_gmm_model_all_w_copy)
        max_weight = nc_ex_gmm_model_all_w_copy[max_weight_index]
        nc_ex_gmm_w_sum += max_weight
        if nc_ex_gmm_w_sum < percentile:
            # Weights we are keeping
            nc_ex_gmm_colour_list[max_weight_index] = colour_in
            nc_ex_gmm_reduced_list.append(nc_ex_gmm_model_all_w_copy[max_weight_index])
            nc_ex_gmm_reduced_labels.append(nc_ex_gmm_model_all_labels[max_weight_index])
            nc_ex_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_labels, \
             'nc_ex_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'Non-cratonic and Extended All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = nc_ex_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    nc_ex_gmm_reduced_list = nc_ex_gmm_reduced_list/sum(nc_ex_gmm_reduced_list)
    bar_plot(nc_ex_gmm_reduced_list, nc_ex_gmm_reduced_list, nc_ex_gmm_reduced_labels, \
             'nc_ex_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'Non-cratonic and Extended Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)

################################
# Sort by weight and identify nth percentile cutoff with reweighting
percentiles = [0.75, 0.8]
colour_in = 'b'
colour_out = 'r'

for percentile in percentiles:
    print 'Cut-off percentile of', (100*percentile)
    banda_gmm_w_sum = 0
    banda_gmm_model_all_w_copy = copy.copy(banda_gmm_model_all_w)
    banda_gmm_colour_list = [colour_out]*len(banda_gmm_model_all_w)
    banda_gmm_reduced_list = [] # list to store models we are keeping
    banda_gmm_reduced_labels = []
    for i in range(len(banda_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(banda_gmm_model_all_w_copy)
        max_weight = banda_gmm_model_all_w_copy[max_weight_index]
        banda_gmm_w_sum += max_weight
        if banda_gmm_w_sum < percentile:
            # Weights we are keeping
            banda_gmm_colour_list[max_weight_index] = colour_in
            banda_gmm_reduced_list.append(banda_gmm_model_all_w_copy[max_weight_index])
            banda_gmm_reduced_labels.append(banda_gmm_model_all_labels[max_weight_index])
            banda_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(banda_gmm_model_all_w, banda_gmm_model_all_w, banda_gmm_model_all_labels, \
             'banda_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'Banda Sea All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = banda_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    banda_gmm_reduced_list = banda_gmm_reduced_list/sum(banda_gmm_reduced_list)
    num_kept_models = len(banda_gmm_reduced_list)
#    print banda_gmm_reduced_list, sum(banda_gmm_reduced_list)
    bar_plot(banda_gmm_reduced_list, banda_gmm_reduced_list, banda_gmm_reduced_labels, \
             'banda_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'Banda Sea Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)

    # Now we want to elimate the lowest weighted models and redistribute within each 
    # region type. This will give a different weighting to the .._reduced_list and
    # aims to preserve the original regional weights
    banda_gmm_orig = copy.copy(banda_gmm_w_list)
    print type(banda_gmm_orig)
    print banda_gmm_orig
    print type(banda_gmm_orig[1])
    num_models = sum(len(item) for item in banda_gmm_orig)
    print num_models
    while num_models > num_kept_models:
        min_weight = 1
        for i in range(len(banda_gmm_orig)):
            if sum(banda_gmm_orig[i]) > 0:
                reg_min = min(k for k in banda_gmm_orig[i] if k > 0)
                reg_min_index = np.where(banda_gmm_orig[i] == reg_min)
                print 'reg_min', reg_min
                print 'reg_min_index', reg_min_index
                if reg_min < min_weight:
                    min_weight = reg_min
                    min_region_index = i
                    j = reg_min_index
                    print 'i,j', i,j
            else:
                pass
        # set minimum element to zero and renormalise region weights
        region_sum = sum(banda_gmm_orig[min_region_index])
        print region_sum
        banda_gmm_orig[min_region_index][j] = 0
        print banda_gmm_orig[min_region_index], sum(banda_gmm_orig[min_region_index])
        if sum(banda_gmm_orig[min_region_index]) > 0:
               banda_gmm_orig[min_region_index] = banda_gmm_orig[min_region_index]/sum(banda_gmm_orig[min_region_index])*region_sum
        print banda_gmm_orig[min_region_index], sum(banda_gmm_orig[min_region_index])
        num_models -= 1
    print banda_gmm_orig#, sum(sum(m for m in banda_gmm_orig))
    banda_gmm_final_weights = np.concatenate(banda_gmm_orig).ravel()
    banda_gmm_final_weights = banda_gmm_final_weights/sum(banda_gmm_final_weights)
    print banda_gmm_final_weights
    print sum(banda_gmm_final_weights)
    bar_plot(banda_gmm_final_weights, banda_gmm_final_weights, banda_gmm_model_all_labels, \
             'banda_gmm_final_weights_%.0fth_percentile.png' % (percentile*100), \
             'Banda Sea Re-normalised Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)
   # sys.exit()

    ###################
    c_gmm_w_sum = 0
    c_gmm_model_all_w_copy = copy.copy(c_gmm_model_all_w)
    c_gmm_colour_list = [colour_out]*len(c_gmm_model_all_w)
    c_gmm_reduced_list = [] # list to store models we are keeping
    c_gmm_reduced_labels = []
    for i in range(len(c_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(c_gmm_model_all_w_copy)
        max_weight = c_gmm_model_all_w_copy[max_weight_index]
        c_gmm_w_sum += max_weight
        if c_gmm_w_sum < percentile:
            # Weights we are keeping
            c_gmm_colour_list[max_weight_index] = colour_in
            c_gmm_reduced_list.append(c_gmm_model_all_w_copy[max_weight_index])
            c_gmm_reduced_labels.append(c_gmm_model_all_labels[max_weight_index])
            c_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(c_gmm_model_all_w, c_gmm_model_all_w, c_gmm_model_all_labels, \
             'c_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'C Sea All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = c_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    c_gmm_reduced_list = c_gmm_reduced_list/sum(c_gmm_reduced_list)
    num_kept_models = len(c_gmm_reduced_list)
#    print c_gmm_reduced_list, sum(c_gmm_reduced_list)
    bar_plot(c_gmm_reduced_list, c_gmm_reduced_list, c_gmm_reduced_labels, \
             'c_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'C Sea Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)

    # Now we want to elimate the lowest weighted models and redistribute within each 
    # region type. This will give a different weighting to the .._reduced_list and
    # aims to preserve the original regional weights
    c_gmm_orig = copy.copy(c_gmm_w_list)
    print type(c_gmm_orig)
    print c_gmm_orig
    print type(c_gmm_orig[1])
    num_models = sum(len(item) for item in c_gmm_orig)
    print num_models
    while num_models > num_kept_models:
        min_weight = 1
        for i in range(len(c_gmm_orig)):
            if sum(c_gmm_orig[i]) > 0:
                reg_min = min(k for k in c_gmm_orig[i] if k > 0)
                reg_min_index = np.where(c_gmm_orig[i] == reg_min)
                print 'reg_min', reg_min
                print 'reg_min_index', reg_min_index
                if reg_min < min_weight:
                    min_weight = reg_min
                    min_region_index = i
                    j = reg_min_index
                    print 'i,j', i,j
            else:
                pass
        # set minimum element to zero and renormalise region weights
        region_sum = sum(c_gmm_orig[min_region_index])
        print region_sum
        c_gmm_orig[min_region_index][j] = 0
        print c_gmm_orig[min_region_index], sum(c_gmm_orig[min_region_index])
        if sum(c_gmm_orig[min_region_index]) > 0:
               c_gmm_orig[min_region_index] = c_gmm_orig[min_region_index]/sum(c_gmm_orig[min_region_index])*region_sum
        print c_gmm_orig[min_region_index], sum(c_gmm_orig[min_region_index])
        num_models -= 1
    print c_gmm_orig#, sum(sum(m for m in c_gmm_orig))
    c_gmm_final_weights = np.concatenate(c_gmm_orig).ravel()
    c_gmm_final_weights = c_gmm_final_weights/sum(c_gmm_final_weights)
    print c_gmm_final_weights
    print sum(c_gmm_final_weights)
    bar_plot(c_gmm_final_weights, c_gmm_final_weights, c_gmm_model_all_labels, \
             'c_gmm_final_weights_%.0fth_percentile.png' % (percentile*100), \
             'Cratonic Re-normalised Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)


    ###########################
    nc_ex_gmm_w_sum = 0
    nc_ex_gmm_model_all_w_copy = copy.copy(nc_ex_gmm_model_all_w)
    nc_ex_gmm_colour_list = [colour_out]*len(nc_ex_gmm_model_all_w)
    nc_ex_gmm_reduced_list = [] # list to store models we are keeping
    nc_ex_gmm_reduced_labels = []
    for i in range(len(nc_ex_gmm_model_all_w_copy)):
        max_weight_index = np.argmax(nc_ex_gmm_model_all_w_copy)
        max_weight = nc_ex_gmm_model_all_w_copy[max_weight_index]
        nc_ex_gmm_w_sum += max_weight
        if nc_ex_gmm_w_sum < percentile:
            # Weights we are keeping
            nc_ex_gmm_colour_list[max_weight_index] = colour_in
            nc_ex_gmm_reduced_list.append(nc_ex_gmm_model_all_w_copy[max_weight_index])
            nc_ex_gmm_reduced_labels.append(nc_ex_gmm_model_all_labels[max_weight_index])
            nc_ex_gmm_model_all_w_copy[max_weight_index] = 0
        else:
            break
    bar_plot(nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_w, nc_ex_gmm_model_all_labels, \
             'nc_ex_gmm_all_weights_%.0fth_percentile.png' % (percentile*100), \
             'Nc_Ex Sea All Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5, \
             colour_list = nc_ex_gmm_colour_list)
    # Re-normalise weights for remaining GMMs
    nc_ex_gmm_reduced_list = nc_ex_gmm_reduced_list/sum(nc_ex_gmm_reduced_list)
    num_kept_models = len(nc_ex_gmm_reduced_list)
#    print nc_ex_gmm_reduced_list, sum(nc_ex_gmm_reduced_list)
    bar_plot(nc_ex_gmm_reduced_list, nc_ex_gmm_reduced_list, nc_ex_gmm_reduced_labels, \
             'nc_ex_gmm_reduced_weights_%.0fth_percentile.png' % (percentile*100), \
             'Nc_Ex Sea Selected Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)

    # Now we want to elimate the lowest weighted models and redistribute within each 
    # region type. This will give a different weighting to the .._reduced_list and
    # aims to preserve the original regional weights
    nc_ex_gmm_orig = copy.copy(nc_ex_gmm_w_list)
    print type(nc_ex_gmm_orig)
    print nc_ex_gmm_orig
    print type(nc_ex_gmm_orig[1])
    num_models = sum(len(item) for item in nc_ex_gmm_orig)
    print num_models
    while num_models > num_kept_models:
        min_weight = 1
        for i in range(len(nc_ex_gmm_orig)):
            if sum(nc_ex_gmm_orig[i]) > 0:
                reg_min = min(k for k in nc_ex_gmm_orig[i] if k > 0)
                reg_min_index = np.where(nc_ex_gmm_orig[i] == reg_min)
                print 'reg_min', reg_min
                print 'reg_min_index', reg_min_index
                if reg_min < min_weight:
                    min_weight = reg_min
                    min_region_index = i
                    j = reg_min_index
                    print 'i,j', i,j
            else:
                pass
        # set minimum element to zero and renormalise region weights
        region_sum = sum(nc_ex_gmm_orig[min_region_index])
        print region_sum
        nc_ex_gmm_orig[min_region_index][j] = 0
        print nc_ex_gmm_orig[min_region_index], sum(nc_ex_gmm_orig[min_region_index])
        if sum(nc_ex_gmm_orig[min_region_index]) > 0:
               nc_ex_gmm_orig[min_region_index] = nc_ex_gmm_orig[min_region_index]/sum(nc_ex_gmm_orig[min_region_index])*region_sum
        print nc_ex_gmm_orig[min_region_index], sum(nc_ex_gmm_orig[min_region_index])
        num_models -= 1
    print nc_ex_gmm_orig#, sum(sum(m for m in nc_ex_gmm_orig))
    nc_ex_gmm_final_weights = np.concatenate(nc_ex_gmm_orig).ravel()
    nc_ex_gmm_final_weights = nc_ex_gmm_final_weights/sum(nc_ex_gmm_final_weights)
    print nc_ex_gmm_final_weights
    print sum(nc_ex_gmm_final_weights)
    bar_plot(nc_ex_gmm_final_weights, nc_ex_gmm_final_weights, nc_ex_gmm_model_all_labels, \
             'nc_ex_gmm_final_weights_%.0fth_percentile.png' % (percentile*100), \
             'Non-cratonic and Extended Re-normalised Model Weights', fig_path = gmm_fig_path, \
             fontsize=10, wide=2.5)
