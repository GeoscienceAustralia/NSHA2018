"""This code combines the smoothed seismicity models applying different b-values and mmax's in different tectonic regions
"""

import os, sys
import numpy as np

from NSHA2018.source_models.logic_trees import logic_tree
from openquake.hazardlib.sourcewriter import write_source_model

def combine_ss_models(filedict, domains_shp, lt, outfile):
    """ Combine smoothed seismicity models based on tectonic region types
    :params filedict:
        dict of form filedict[trt] = filename specifying input file for that region
    :params domains_shp:
        shapefile defining tectonic domain regions
    :params lt:
        LogicTree object containing relevant values and weights for Mmax
    :params outfile:
        output nrml formatted file
    """







if __name__ == "__main__":
    filedict = {'Cratonic': 'source_model_Australia_Adaptive_K3_b0.819.csv.xml',
                'Non-Cratonic': 'source_model_Australia_Adaptive_K3_b1.208.csv.xml',
                'Extended': 'source_model_Australia_Adaptive_K3_b0.835.csv.xml'}
    domains_shp = '../zones/2012_mw_ge_4.0/NSHA13_Background/shapefiles/NSHA13_BACKGROUND_NSHA18_MFD.shp'
    outfile = 'source_model_Australia_Adaptive_K3_merged.xml'
    lt  = logic_tree.LogicTree('../../shared/seismic_source_model_weights_rounded_p0.4.csv')
    combine_ss_models(filedict, domains_shp, lt, outfile)
