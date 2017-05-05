"""Plot a series of hazard curves from OpenQuake output

Jonathan Griffin
Geoscience Australia, January 2017
"""

import matplotlib as mpl
mpl.use('Agg')
import os, sys
import matplotlib.pyplot as plt

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
                if float(row[0]) == location[0] and float(row[1]) == location[1]:
#                    annual_probs.append([row[2:]])
                    plt.semilogy(PGA_values, row[2:])
                    print row[2:]
        if labels is not None:
            plt.legend(labels, loc = 'upper right')
        plt.xlabel('PGA(g)')
        plt.ylabel('POE in 50 years')
        plt.xlim(0,2.5)
        plt.ylim(10E-7,10E-2)
        plt.savefig('hazard_curves_%.6f_%.6f.png' % (location[0], location[1]))
        plt.clf()


if __name__ == "__main__":
    file_list = []
#    pathname = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_runs/Cadell_testing/Cadell_cluster'
 #   pathname = '/short/w84/NSHA18/PSHA_modelling/Cadell_testing/Cadell_cluster/'
    pathname = '/short/w84/NSHA18/PSHA_modelling'
#    folder_list = ['results_cadell_cluster_32000',
#                   'results_cadell_cluster_34500',
#                   'results_cadell_cluster_background',
#                   'results_cadell_cluster_start']

    folder_list = ['Adelaide_zone_as_pts/20170503_154940_Adelaide_zone_as_pts_jdg547/results',
                   'Adelaide_zone/20170503_135847_Adelaide_zone_jdg547/results']
#    locations = [[144.77539, -35.71946]]
    
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
    labels = ['Zone as pts', 'Zone']
    plot_hazard_curves(file_list, locations, labels = labels)

