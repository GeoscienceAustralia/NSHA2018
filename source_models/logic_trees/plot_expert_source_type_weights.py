"""Plot the range of response for each experts assessment of the weight for each source model type
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

# Load in target responses for seismic source mode,l
target_path = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHA_18/Expert_Elicitation/Target_questions'

SS_weights = []
BG_weights = []
Reg_weights = []
Seismo_weights = []
SS_flt_weights = []

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
        S2 = ids[35:40]

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
    S2_sum = sum(data[35:40,1])
    print data[35]
    SS_weights.append(data[35])
    BG_weights.append(data[36])
    Reg_weights.append(data[37])
    Seismo_weights.append(data[38])
    SS_flt_weights.append(data[39])

SS_weights = np.array(SS_weights)
BG_weights = np.array(BG_weights)
Reg_weights = np.array(Reg_weights)
Seismo_weights = np.array(Seismo_weights)
SS_flt_weights = np.array(SS_flt_weights)
print SS_weights
print BG_weights
print Reg_weights
print Seismo_weights
print SS_flt_weights
weight_list = [SS_weights, BG_weights, Reg_weights, Seismo_weights,
               SS_flt_weights]
fig = plt.figure(figsize = (3,10))
title_list = ['Smoothed Seismicity', 'Background', 'Regional',
              'Seismotectonic', 'Smoothed Seismicity with Faults']
label_list = ['a)', 'b)', 'c)', 'd)', 'e)']
for i,mod_type_w in enumerate(weight_list):
    plt.subplot(5,1,i+1)
    for j in range(num_ssc_experts):
        plt.scatter(mod_type_w[j,1], [j+1], marker='x', c='k')
        plt.plot([mod_type_w[j,0], mod_type_w[j,2]], [j+1, j+1], c='k')
    plt.yticks(np.arange(1, num_ssc_experts+1, 1), fontsize=7.5)
    plt.xticks(np.arange(0.0, 1.2, 0.2), fontsize=7.5)
    plt.xlim([-0.05,1.0])
    plt.ylim([0, 17])
    plt.xlabel('Weight', fontsize=10)
    plt.ylabel('Expert ID', fontsize=10)
    plt.text(0, 15.5, label_list[i], fontsize=10)
#    plt.title(title_list[i], fontsize=10)

plt.tight_layout()
plt.savefig('Model_type_expert_weights.png', dpi=300)
plt.savefig('Model_type_expert_weights.eps', dpi=300)
