# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 14:39:20 2019

@author: u93322
"""
#import IPython
#IPython.embed()

import glob
import sys, os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from palettable.colorbrewer.qualitative import Set1_5

def main():

#TODO output a file of scenarios:
#     :name, location, mag max, loc max, poe
#TODO make sure paths work on both unix and windows - may be easier to run on NCI when finalised?
    
    place_name = "Newcastle" # TODO - possibility for being a list of places in file
    poe = 0.005
    df_mag = open_rlz_file(place_name, poe, "mag")
    top_mag = select_top_mags(df_mag)
    print("The highest probability magnitude for poe: %s, is %s: " % (poe, top_mag))
    print("looking up modal location...")
    df_loc = open_rlz_file(place_name, poe, "mll")
    max_loc = select_top_loc(df_loc, top_mag)
    details = append_to_scenrio_list(max_loc, place_name)
    print(details)

def open_rlz_file(place_name, poe, deagg_type):

    file_path = Path("%s/%s" % (place_name, deagg_type))
    for file in os.listdir(str(file_path)):
        if file.endswith(".csv"):
            path_to_file = file_path / str(file)
            with open(str(path_to_file)) as f:
                first_line = f.readline()
                headers = first_line.split(', ')
                poe_file = headers[6].split('=')[1]
                poe_file = np.float(re.sub("'", '', poe_file))
                if poe_file == poe:
                    infile = str(path_to_file)
                    df = make_dataframes(infile)
    return df


def make_dataframes(infile):
    df = pd.read_csv(infile, header=1)
    # replace blank spaces with "_"
    df.columns = [column.replace(" ", "_") for column in df.columns]
    
    #check for missing data
    if df.isnull().sum().any() != 0:
        print(df.isnull().sum())
        sys.exit("Missing values in data file... exiting..")
    return df

            
def select_top_mags(df):
    max_poe = df.nlargest(1, ['poe'])
    mag_max = np.float(max_poe.mag)
    return mag_max


def select_top_loc(df, max_mag, plot=True):
    #fix wrongly labeled columns from OQ
    if df.columns[0] == "mag":
        df.rename(columns = {"mag": "lon", 
                             "lon":"lat",
                             "lat":"mag"}, 
                             inplace = True)
    locs = df.loc[df['mag'] == max_mag]
    if plot == True:
        locs.plot.line(y='poe')
        plt.show()
    max_loc = locs.loc[locs['poe'].idxmax()]
    return max_loc


def append_to_scenrio_list(df, place_name):
    file_line = [place_name, df['mag'], df['lat'], df['lon'], df['poe']]
    return file_line


if __name__ == "__main__":
    main()
