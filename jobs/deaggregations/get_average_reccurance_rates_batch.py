import os, sys
import csv
import platform
import numpy as np
import pandas as pd
from pathlib import Path
from glob import glob
import shapefile
from shapely.geometry import Point, Polygon

"""
Script used to calcualte recurrance intervals and their corresponding 
uncertainty ranges.  Used specifically for creation of the Altas of
Earthquake Scenarios.
"""

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

# TODO calculate mean/stdev from mfds in log space.  

def main():

    # Make sure keys are EXACTLY the same as the NSHA folder names
    source_model_dict = {
                        "AUS6": 0.204,
                        "DIMAUS": 0.210,
                        "NSHA13": 0.283,
                        "ARUP": 0.052,
                        "ARUP_Background": 0.052,
                        "Leonard2008": 0.053,
                        "Domains_multi_mc": 0.092,
                        "NSHA13_Background": 0.036,
                        "SinMcC2016": 0.019
    }

    # from file or hard coded
    # TODO put as a choice between a looping file or single input arguments
    reccurance_list = []

    with open("Earthquake_Scenarios2.csv") as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            name = row[0]
            lon, lat = np.float(row[3]), np.float(row[2])
            mag = np.float(row[1])

            x_vals, weights = look_up_recurrance(lon, lat, mag, source_model_dict)            
            mean, stdev = weighted_avg_and_std(x_vals, weights)

            #print(mean)
            #print(stdev)
            ub = mean - (2*stdev)
            lb = mean + (2*stdev)

            ub_exp = np.exp(ub)
            lb_exp = np.exp(lb)

            lb_yrs = (1 / lb_exp)
            ub_yrs = (1 / ub_exp)
            
            mean_exp = np.exp(mean)
            mean_yrs = np.round(1 / mean_exp)
            
            print("For Scenario %s: " % name)
            print("mean recurrance = %s in the range %s to %s" % (mean_yrs, np.round(lb_yrs), np.round(ub_yrs)))
            
            reccurance_list.append([name, mean_yrs, np.round(lb_yrs), np.round(ub_yrs)])


    np.savetxt(u"Earthquake_Scenarios2_rates.csv", np.array(reccurance_list), fmt=u'%10.20s', delimiter=u",")
    #print(reccurance_list)    


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    vals_log = np.log(values)
    wtaverage = np.average(vals_log, weights=weights)
    variance = np.average((vals_log-wtaverage)**2, weights=weights)
    return (wtaverage, np.sqrt(variance))


def look_up_recurrance(lon, lat, mag, source_model_dict):
    
    weights = []
    x_vals = []

    for key in source_model_dict.keys():
        # relative path for source model 
  
        NSHA_SM_path = Path("../../source_models/zones/2018_mw/")
        full_path = NSHA_SM_path / str(key) / "shapefiles"

        for file in os.listdir(str(full_path)):
            if file.endswith(".shp"):
                sf = shapefile.Reader(str(full_path / str(file)))
        shapes = sf.shapes()
        records = sf.records()

        for shape, record in zip(shapes, records):
            polygon = Polygon(shape.points)
            point = Point(lon, lat)
            if point.within(polygon):
                zone_code = record.CODE

        #print("key = %s" % key)
        #print("zone code = %s" % zone_code )

        rates_file = str(zone_code) + "_rates.csv"
        rates_path_file = NSHA_SM_path / str(key) / "mfd" / str(zone_code) / rates_file
        
        df = pd.read_csv(rates_path_file)
        # get value closest in magnitude as these mdf files arem't exactly the same
        line = df.iloc[(df['MAG']-mag).abs().argsort()[:1]]

        mfd_fit_area_norm = line['MFD_FIT_AREA_NORM'].item()

        weights.append(source_model_dict[key])
        x_vals.append(mfd_fit_area_norm)

    return x_vals, weights

if __name__ == "__main__":
    main()
