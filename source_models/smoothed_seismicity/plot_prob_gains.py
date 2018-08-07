"""Script for making plots of probability gain statistics
"""

import os,sys
from os.path import join
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

folder = 'llh_results'
infiles_K4_3_5_L65 = ['Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin3.5_learning1965_1974'
xaxis = 'target' # 'learning' or 'target'

infiles_K3_4_L65 = ['Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin4.0_learning1965_1974'
xaxis = 'target' # 'learning' or 'target

infiles_K4_4_L65 = ['Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin4.0_learning1965_1974'
xaxis = 'target' # 'learning' or 'target 

infiles_K5_4_L65 = ['Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin4.0_learning1965_1974'
xaxis = 'target' # 'learning' or 'target 

infiles_K3_4 = ['Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin4.0_target2005_2010'
xaxis = 'learning' # 'learning' or 'target'

infiles_K4_4 = ['Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin4.0_target2005_2010'
xaxis = 'learning'

infiles_K5_4 = ['Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin4.0_target2005_2010'
xaxis = 'learning'


#infiles = ['Australia_Adaptive_K4_res0.50_b1.000_mmin3.5_learning1965_1974_target1975_1984_llh.csv',
#           'Australia_Adaptive_K4_res0.50_b1.000_mmin3.5_learning1965_1974_target1985_1994_llh.csv',
#           'Australia_Adaptive_K4_res0.50_b1.000_mmin3.5_learning1965_1974_target1995_2004_llh.csv',
#           'Australia_Adaptive_K4_res0.50_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv']
#figure_comment = 'K4_res0.50_b1.000_mmin3.5_learning1965_1974'
#xaxis = 'target' # 'learning' or 'target'

infiles_K3_3_5 = ['Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin3.5_target2005_2010'
xaxis = 'learning'

infiles_K4_3_5 = ['Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin3.5_target2005_2010'
xaxis = 'learning' # 'learning' or 'target'

infiles_K5_3_5_L65 = ['Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin3.5_learning1965_1974'
xaxis = 'target' # 'learning' or 'target'

infiles_K3_3_5_L65 = ['Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_1984_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target1985_1994_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target1995_2004_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin3.5_learning1965_1974'
xaxis = 'target' # 'learning' or 'target'

infiles_K5_3_5 = ['Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1975_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1985_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1995_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin3.5_target2005_2010'
xaxis = 'learning' # 'learning' or 'target'

infiles_K3_4_X = ['Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'

infiles_K3_3_5_X = ['Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'

infiles_K4_4_X = ['Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'

infiles_K4_3_5_X = ['Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'

infiles_K5_4_X = ['Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'

infiles_K5_3_5_X = ['Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target1975_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1984_target1985_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1994_target1995_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'target' # 'learning' or 'target'



infiles_K3_4_Y = ['Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'

infiles_K3_3_5_Y = ['Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K3_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K3_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'

infiles_K4_4_Y = ['Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'

infiles_K4_3_5_Y = ['Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K4_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K4_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'

infiles_K5_4_Y = ['Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin4.0_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin4.0_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'

infiles_K5_3_5_Y = ['Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1974_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1984_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_1994_target2005_2010_llh.csv',
           'Australia_Adaptive_K5_res0.10_b1.000_mmin3.5_learning1965_2004_target2005_2010_llh.csv']
figure_comment = 'K5_res0.10_b1.000_mmin3.5_learning1965_1x'
xaxis = 'learning' # 'learning' or 'target'


#infiles_list = [infiles_K3_3_5, infiles_K4_3_5, infiles_K5_3_5, infiles_K3_4, infiles_K4_4, infiles_K5_4]
#labels = ['K=3, MMin=3.5', 'K=4, MMin=3.5', 'K=5, MMin=3.5','K=3, MMin=4.0', 'K=4, MMin=4.0', 'K=5, MMin=4.0']
#figure_comment='K3_4_Mmin_3p5_4_target_2005_2010'
#title = 'Probability gain for target period 2005-2010'
#title = 'Kagan I0 for target period 2005-2010'
#title = 'Log likelihood of smoothed model for target period 2005-2010'
#title = 'Log likelihood of uniform model for target period 2005-2010'
#xaxis='learning'
#loc=2 # upper left
#loc=1 #upper right
#loc=4 #lower right
#loc=6 #centre left
#loc=7 #centre right

#infiles_list = [infiles_K3_3_5_L65, infiles_K4_3_5_L65, infiles_K5_3_5_L65, infiles_K3_4_L65, infiles_K4_4_L65, infiles_K5_4_L65]
#labels = ['K=3, MMin=3.5', 'K=4, MMin=3.5', 'K=5, MMin=3.5','K=3, MMin=4.0', 'K=4, MMin=4.0', 'K=5, MMin=4.0']
#figure_comment='K3_4_Mmin_3p5_4_learning_1967_1974'
#title = 'Probability gain for learning period 1965-1974'
#title = 'Kagan I1 for learning period 1965-1974'
#loc=2
#title = 'Log likelihood of uniform model for learning period 1965-1974'
#title = 'Log likelihood of smoothed model for learning period 1965-1974'
#xaxis='target'
#loc=1
#loc=4
#loc=6

#infiles_list = [infiles_K3_3_5_X, infiles_K4_3_5_X, infiles_K5_3_5_X, infiles_K3_4_X, infiles_K4_4_X, infiles_K5_4_X]
#labels = ['K=3, MMin=3.5', 'K=4, MMin=3.5', 'K=5, MMin=3.5','K=3, MMin=4.0', 'K=4, MMin=4.0', 'K=5, MMin=4.0']
#figure_comment='K3_4_Mmin_3p5_4_learning_1965_X'
#title = 'Probability gain for learning period 1965-X \n against remainder of catalogue'

#title = 'Kagan I1 for learning period 1965-1974'
#loc=2
#title = 'Log likelihood of uniform model for learning period 1965-1974'
#title = 'Log likelihood of smoothed model for learning period 1965-1974'
#xaxis='target'
#loc=1
#loc=4
#loc=2

infiles_list = [infiles_K3_3_5_Y, infiles_K4_3_5_Y, infiles_K5_3_5_Y, infiles_K3_4_Y, infiles_K4_4_Y, infiles_K5_4_Y]
labels = ['K=3, MMin=3.5', 'K=4, MMin=3.5', 'K=5, MMin=3.5','K=3, MMin=4.0', 'K=4, MMin=4.0', 'K=5, MMin=4.0']
#infiles_list = [infiles_K3_3_5_Y, infiles_K3_4_Y]
#labels = ['MMin=3.5','MMin=4.0']
figure_comment='K3_4_Mmin_3p5_4_learning_1965_Y'
title = 'Probability gain for learning period 1965-X \n and target period 2005-2010'
#title = 'Kagan I1 for learning period 1965-1974'
loc=2
#title = 'Log likelihood of uniform model for learning period 1965-1974'
#title = 'Log likelihood of smoothed model for learning \n period 1965-X and target period 2005-2010'
xaxis='learning'
#loc=1
#loc=4
#loc=2

markers = ['o', '^', 'p', 'h', 'D', 'x']
colours = [(1.0, 0., 0.), (1.0, 0.7, 0.), (1.0, 1.0, 0.),'c', 'g', (0., 1.0, 0.05)]

for i,infiles in enumerate(infiles_list):
    data = []
    #labels = []
    target_dates = []
    learning_dates = []
    target_width = []
    learning_width = []
    for filename in infiles:
	try:
	    learning_period_start = filename.split('_')[5].lstrip('learning')
	    learning_period_end = filename.split('_')[6]
	    print (float(learning_period_end) - float(learning_period_start))
	    print (float(learning_period_end) - (float(learning_period_end) - float(learning_period_start)))
            learning_width.append((float(learning_period_end) - float(learning_period_start))/2)
	    learning_dates.append(float(learning_period_end) - ((float(learning_period_end) - \
								     float(learning_period_start))/2))
	    learning_period = learning_period_start + '-' + learning_period_end

	    target_period_start = filename.split('_')[7].lstrip('target')
	    target_period_end = filename.split('_')[8]
            target_width.append((float(target_period_end) - float(target_period_start))/2)
	    target_dates.append(float(target_period_end) - (float(target_period_end) - \
								float(target_period_start)))
	    target_period = target_period_start + '-' + target_period_end
	except ValueError:
	    learning_period_start = filename.split('_')[6].lstrip('learning')
	    learning_period_end = filename.split('_')[7]
	    print (float(learning_period_end) - float(learning_period_start))
	    print (float(learning_period_end) - (float(learning_period_end) - float(learning_period_start)))
            learning_width.append(((float(learning_period_end) - float(learning_period_start))+1)/2)
	    learning_dates.append(float(learning_period_end) - ((float(learning_period_end) - \
								     float(learning_period_start))/2))
	    learning_period = learning_period_start + '-' + learning_period_end

	    target_period_start = filename.split('_')[8].lstrip('target')
	    target_period_end = filename.split('_')[9]
            width = ((float(target_period_end) - float(target_period_start))+1)/2
            target_width.append(width)
	    target_dates.append(float(target_period_end) - width)
	    target_period = target_period_start + '-' + target_period_end

	#labels.append(target_period)
	filepath = os.path.join(folder, filename)
	tmp_data = np.genfromtxt(filepath, delimiter=',')
	data.append(tmp_data)
    data = np.array(data)
    print data
    print data[:,4]
    if xaxis=='learning':
	xs = learning_dates
        xerr = learning_width
    elif xaxis=='target':
	xs = target_dates
        xerr = target_width
    print xs
    print xerr

    #plt.plot(xs, -1*data[:,0])
   # plt.scatter(xs, -1*data[:,0],marker=markers[i], c= colours[i], s=60)
   # plt.errorbar(xs, -1*data[:,0], xerr=xerr, linestyle='None', c=colours[i], capsize=0)
    #plt.savefig(join(folder,('model_log_likelihoods_%s.png' % figure_comment)))
#    plt.clf()
    #plt.plot(xs, data[:,1])
    #plt.scatter(xs, data[:,1],marker=markers[i], c= colours[i], s=60)
    #plt.errorbar(xs, data[:,1], xerr=xerr, linestyle='None', c=colours[i], capsize=0)
    #plt.savefig(join(folder,('kagan_IO_%s.png' % figure_comment)))
    #plt.clf()
    #plt.plot(xs, data[:,2])
    #plt.scatter(xs, data[:,2],marker=markers[i], c= colours[i], s=60)
    #plt.errorbar(xs, data[:,2], xerr=xerr, linestyle='None', c=colours[i], capsize=0)
    #plt.savefig(join(folder,('kagan_I1_%s.png' % figure_comment)))
    #plt.clf()
    #plt.plot(xs, data[:,3])
    #plt.scatter(xs, data[:,3],marker=markers[i], c= colours[i], s=60)
    #plt.errorbar(xs, data[:,3], xerr=xerr, linestyle='None', c=colours[i], capsize=0)
    #plt.savefig(join(folder,('uniform_log_likelihoods_%s.png' % figure_comment)))
    #plt.clf()
    # probability gains
    plt.scatter(xs, data[:,4],marker=markers[i], c= colours[i], s=60)
    plt.errorbar(xs, data[:,4], xerr=xerr, linestyle='None', c=colours[i], capsize=0)
   # plt.scatter(xs, data[:,4],marker=markers[i], c= colours[i])
plt.ylabel('Probability Gain')
#plt.ylabel('Information Score')
#plt.ylabel('Log Likelihood')
if xaxis=='learning':
    plt.xlabel('Learning Period')
elif xaxis=='target':
    plt.xlabel('Target Period')
plt.title(title)

if figure_comment=='K3_4_Mmin_3p5_4_learning_1965_Y':
    lgd = plt.legend(labels, loc=9, scatterpoints=1, bbox_to_anchor=(0.5,-0.1), ncol=3)
    #plt.savefig(join(folder,('probability_gains_%s.png' % figure_comment)), resolution=900,
    #        bbox_extra_artists=(lgd,), bbox_inches='tight' )
    plt.savefig(join(folder,('probability_gains_%s.png' % figure_comment)), resolution=900,
                bbox_extra_artists=(lgd,), bbox_inches='tight' ) 
else:
    lgd = plt.legend(labels, loc=loc, scatterpoints=1)#, ncol=2)
    plt.savefig(join(folder,('probability_gains_%s.png' % figure_comment)), resolution=900,
                bbox_inches='tight' )
#plt.savefig(join(folder,('kagan_I1_%s.png' % figure_comment)), resolution=900)
#plt.savefig(join(folder,('model_log_likelihoods_%s.png' % figure_comment)), resolution=900)
