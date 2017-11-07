"""Plot a series of hazard curves from OpenQuake output

Jonathan Griffin
Geoscience Australia, January 2017
"""

import matplotlib as mpl
mpl.use('Agg')
import os, sys
import matplotlib.pyplot as plt
import numpy as np

def plot_hazard_curves(file_list, location_list, labels = None):
    """Plot hazard curves on a single figure by looping through
    all provided files and locations
    :param file_list:
        list of filepaths containing hazard curves in OpenQuake
        csv output format
    :param location_list:
        list of pairs of [[latitude, longitude]]. Must exist in the
        hazard curve files
    """

    # Read the data
    for location in location_list:
        print location
        for filepath in file_list:
            #print filepath
            f_in = open(filepath, 'r')
            header1 = f_in.readline() # ignore
            header2 = f_in.readline().split(',')
            PGA_values = []
            PGA_names = header2[2:]
            for PGA in PGA_names:
                value = float(PGA.lstrip('PGA-'))
                PGA_values.append(value)
            annual_probs = []
            for line in f_in.readlines(): # Now read through data to find points of interest
                row = line.split(',')
            #    print row
#                for location in location_list:
#                print location
#                print float(row[0]), float(row[1])
                if np.allclose(float(row[0]), location[0], rtol=0.001) and np.allclose(float(row[1]), location[1], rtol=0.001):
#                    annual_probs.append([row[2:]])
                    plt.semilogy(PGA_values, row[2:])
                    print row[2:]
        if labels is not None:
            plt.legend(labels, loc = 'upper right')
        plt.xlabel('PGA(g)')
        plt.ylabel('POE in 50 years')
        plt.xlim(0,0.5)
        plt.ylim(10E-5,10E-0)
        print filepath
        figname = filepath[:-4] + '_hazard_curves_%.6f_%.6f.png' % (location[0], location[1])
        plt.savefig(figname)
        plt.clf()


if __name__ == "__main__":
    file_list = []
#    pathname = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_runs/Cadell_testing/Cadell_cluster'
 #   pathname = '/short/w84/NSHA18/PSHA_modelling/Cadell_testing/Cadell_cluster/'
    pathname = '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed'
#    folder_list = ['results_cadell_cluster_32000',
#                   'results_cadell_cluster_34500',
#                   'results_cadell_cluster_background',
#                   'results_cadell_cluster_start']

#    folder_list = ['Adelaide_zone_as_pts/20170503_154940_Adelaide_zone_as_pts_jdg547/results',
#                   'Adelaide_zone/20170503_135847_Adelaide_zone_jdg547/results']#    folder_list = ['NSHA_Production_May17/smoothed/GA_adaptive_smoothing/best_bval/hazcurve_20170512_093334__jdg547/results_localities']
#                   'Adelaide_zone/20170503_135847_Adelaide_zone_jdg547/results']
#    locations = [[144.77539, -35.71946]]
#    folder_list = ['NFSM_Collapsed_testing/NFSM_collapsed_additive_w1/hazcurve_20170515_124344_NFSM_collapsed_additive_w1_jdg547/results_localities',
#                   'NFSM_AUS6_Collapsed_testing/NFSM_AUS6_collapsed_additive_w1/hazcurve_20170515_125304_NFSM_AUS6_collapsed_additive_w1_jdg547/results_localities',
#                   'NFSM_DIMAUS_Collapsed_testing/NFSM_DIMAUS_collapsed_additive_w1/hazcurve_20170515_125940_NFSM_DIMAUS_collapsed_additive_w1_jdg547/results_localities']
    
    folder_list = ['/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed/GA_adaptive_smoothing_collapsed/b_mmax_uncert/GA_adaptive_smoothing_collapsed_K4_mmin3p0/hazcurve_20170523_061457_GA_adaptive_smoothing_collapsed_K4_mmin3p0_jdg547/results_localities',
                   '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed/GA_adaptive_smoothing_collapsed/b_mmax_uncert/GA_adaptive_smoothing_collapsed_K4_mmin3p5/hazcurve_20170523_061434_GA_adaptive_smoothing_collapsed_K4_mmin3p5_jdg547/results_localities',
                   '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/smoothed/GA_adaptive_smoothing_collapsed/b_mmax_uncert/GA_adaptive_smoothing_collapsed_K4_mmin4p0/hazcurve_20170523_061450_GA_adaptive_smoothing_collapsed_K4_mmin4p0_jdg547/results_localities',
                   '/short/w84/NSHA18/PSHA_modelling/NSHA_Production_May17/complete_model/no_smoothed/hazcurve_20170521_113733_complete_model_tia547/results_localities']

    for folder in folder_list:
        file_list.append(os.path.join(pathname, folder, 'hazard_curve-mean_1.csv'))
    # get all hazard point locations
    locations =[]    
    f_in = open(file_list[0], 'r')
    header1 = f_in.readline()
    header2 = f_in.readline()
    for line in f_in.readlines():
        row = line.split(',')
        locations.append([float(row[0]), float(row[1])])
    print locations
#    labels = ['Present rate', 'Present + 2500 years rate', 'Background rate', 'Start of cluster rate']
#    labels = ['NSHA13_seismo', 'AUS6_seismo', 'DIMAUS_seismo']
    labels = ['AdapK4Mmin3.0', 'AdapK4Mmin3.5', 'AdapK4Mmin4.0','FullModel_no_ss' ]
#    labels = ['Adaptive_collapsed', 'Fixed_collapsed'] 
    plot_hazard_curves(file_list, locations, labels = labels)

