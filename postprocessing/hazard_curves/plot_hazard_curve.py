"""Plot a series of hazard curves from OpenQuake output

Jonathan Griffin
Geoscience Australia, January 2017
"""

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
    for filepath in file_list:
        print filepath
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
            for location in locations:
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
    plt.savefig('hazard_curves.png')
        


if __name__ == "__main__":
    file_list = []
    pathname = '/nas/gemd/ehp/georisk_earthquake/hazard/NSHM_18/PSHA_modelling/oq_runs/Cadell_testing/Cadell_cluster'
    folder_list = ['results_cadell_cluster', 
                   'results_cadell_cluster_background',
                   'results_cadell_cluster_start']
    locations = [[144.77539, -35.71946]]
    for folder in folder_list:
        file_list.append(os.path.join(pathname, folder, 'hazard_curve-mean_1.csv'))
    labels = ['Present rate', 'Background rate', 'Start of cluster rate']
    plot_hazard_curves(file_list, locations, labels = labels)

