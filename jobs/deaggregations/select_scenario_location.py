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
from os.path import join

import numpy as np
import pandas as pd

def main():

#TODO Add function to calculate the distance of the scenario location to the population centre
#TODO Decide if normalsed mag needed.  

    # open pre canned paths of run deaggregations (set is final run paths)
    f = open("Job_list_20200324_145200.txt", 'r')
    #scenario_list = ["city", "mag", "lat", "lon", "poe", "poe_mag"]
    scenario_list1 = []
    scenario_list2 = []

    for line in f:
    #print(line)
        print("**********")
        line = re.sub('\n', '', line).strip() # remove \n characters
        print(line)
        city = line.split("/")[6]
        full_line = join(line, "results")
        poe = float(0.005)

        df = open_rlz_file(city, poe, full_line)

        # stack poes in each lon/lat bin
        mags, lonlats_df = stack_mags_lonlat(df)

        # select highest poe from lonlat df
        lon, lat = select_top_loc(lonlats_df)

        # make dataframe for only selected lon and lats 
        max_loc_df = df.loc[(df['lon'] == float(lon)) & (df['lat'] == float(lat))]

        # look up highest poe mags for this location
        top_mag1 = select_top_mags(max_loc_df, number=1, mag_prev="NaN")
        top_mag2 = select_top_mags(max_loc_df, number=2, mag_prev="NaN")

        sum_mag1 = sum_locs(df, top_mag1)
        sum_mag2 = sum_locs(df, top_mag2)

        #details1 = append_to_scenrio_list(lon, lat, city, sum_mag1)        
        #details2 = append_to_scenrio_list(lon lat, city, sum_mag2)

        details1 = [city, top_mag1, float(lat), float(lon), sum_mag1]
        details2 = [city, top_mag2, float(lat), float(lon), sum_mag2]

        scenario_list1.append(details1)
        scenario_array1 = np.array(scenario_list1)

        scenario_list2.append(details2)
        scenario_array2 = np.array(scenario_list2)

        np.savetxt(u"Earthquake_Scenarios1.csv", scenario_array1, fmt=u'%10.20s', delimiter=u",")
        np.savetxt(u"Earthquake_Scenarios2.csv", scenario_array2, fmt=u'%10.20s', delimiter=u",")


def open_rlz_file(place_name, poe, file_path):

    for file in os.listdir(str(file_path)):
        if file.startswith("rlz"):
            path_to_file = os.path.join(file_path, str(file))

            with open(str(path_to_file)) as f:
                first_line = f.readline()
                headers = first_line.split(', ')
                poe_file = headers[9].split('=')[1]
                poe_file = np.float(re.sub("'", '', poe_file))

                if poe_file == poe:
                    infile = str(path_to_file)
                    df = make_dataframes(infile)

                    #fix wrongly labeled columns from OQ
                    if df.columns[0] == "mag":
                        df.rename(columns = {"mag": "lon", 
                                "lon":"lat",
                                "lat":"mag"}, 
                                inplace = True)
                    return df


def make_dataframes(infile):
    import pandas as pd
    df = pd.read_csv(infile, header=1)
    # replace blank spaces with "_"
    df.columns = [column.replace(" ", "_") for column in df.columns]
    
    #check for missing data
    if df.isnull().sum().any() != 0:
        print(df.isnull().sum())
        sys.exit("Missing values in data file... exiting..")
    return df


def select_top_mags(df, number=1, mag_prev="NaN"):
    if number == 1:
        max_poe = df.nlargest(1, ['poe'])
        mag_max = np.float(max_poe.mag)
        print(mag_max)
        return mag_max
    if number >= 2:
        i = number
        #print("number is currently %s" % i)
        max_poe = df.nlargest(i, ['poe']).iloc[[i-1]]
        mag_max = np.float(max_poe.mag)
        mag_prev = df.nlargest(1, ['poe']).mag
        if abs(float(mag_max) - float(mag_prev)) < 0.5:
            number = i + 1
            return select_top_mags(df, number, mag_max)
        else:
            if mag_max:
                return mag_max


def stack_mags_lonlat(df):
    #get unique magnitudes
    mags = np.unique(df.mag)
    p_norm = df.poe / sum(df.poe)

    #get stacked sum of z's for a given lat/lon
    poe_sum = []
    for i in range(0, len(p_norm), len(mags)):
        poe_sum.append(sum(p_norm[i:i+len(mags)]))
        
    lonlats_df = df[['lon', 'lat']].drop_duplicates().reset_index()
    lonlats_df['poe_sum'] = poe_sum

    return mags, lonlats_df


def select_top_loc(df):
    max_lonlat = df.nlargest(1, ['poe_sum']) 
    #max_loc = locs.loc[locs['poe'].idxmax()]
    return max_lonlat.lon, max_lonlat.lat


def append_to_scenrio_list(df, place_name, sum_mag):
    file_line = [place_name, df['mag'], df['lat'], df['lon'], df['poe'], sum_mag]
    return file_line


def sum_locs(df, mag):
    locs = df.loc[df['mag'] == mag]
    sum_mag = locs['poe'].sum()
    return sum_mag




if __name__ == "__main__":
    main()
