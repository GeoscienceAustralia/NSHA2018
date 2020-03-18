import os, sys
import platform
import numpy as np
import pandas as pd
from pathlib import Path
from glob import glob
import shapefile
from shapely.geometry import Point, Polygon

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
    lon, lat = 145.1, -37.9
    mag = 6.25

    x_vals, weights = look_up_recurrance(lon, lat, mag, source_model_dict)            
    mean, stdev = weighted_avg_and_std(x_vals, weights)
    print(mean, stdev)
    print(1/mean, stdev)


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    wtaverage = np.average(values, weights=weights)
    variance = np.average((values-wtaverage)**2, weights=weights)
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

        rates_file = str(zone_code) + "_rates.csv"
        rates_path_file = NSHA_SM_path / str(key) / "mfd" / str(zone_code) / rates_file
        
        df = pd.read_csv(rates_path_file)
        line = df.loc[(df['MAG'] == mag)]
        mfd_fit_area_norm = line['MFD_FIT_AREA_NORM'].item()

        weights.append(source_model_dict[key])
        x_vals.append(mfd_fit_area_norm)

    return x_vals, weights

if __name__ == "__main__":
    main()





