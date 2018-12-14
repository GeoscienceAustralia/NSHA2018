# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:01:44 2018

@author: u56903
"""

import shapefile
from sys import argv
from tools.nsha_tools import get_field_data
from tools.source_shapefile_builder import get_preferred_catalogue, \
                                           get_completeness_model, get_aus_shmax_vectors, \
                                           get_rate_adjust_factor, build_source_shape, \
                                           get_ul_seismo_depths, get_neotectonic_domain_params, \
                                           aggregate_intraslab_sources

shpPath = argv[1]

# read shapefile
print 'Reading source shapefile...'
sf = shapefile.Reader(shpPath)
shapes = sf.shapes()
fields = sf.fields

'''
newFields = []
fieldType = []
fieldSize = []
fieldDecimal = []
for f in fields[1:]:
    newFields.append(f[0])
    fieldType.append(f[1])
    fieldSize.append(f[2])
    fieldDecimal.append(f[3])
'''

newFields = ['SRC_NAME', 'CODE', 'SRC_TYPE', 'CLASS', 'SRC_WEIGHT', 'RTE_ADJ_F',
               'DEP_BEST', 'DEP_UPPER', 'DEP_LOWER', 'USD', 'LSD', 'OW_LSD',
               'MIN_MAG', 'MIN_RMAG', 'MMAX_BEST', 'MMAX_LOWER', 'MMAX_UPPER',
               'N0_BEST', 'N0_LOWER', 'N0_UPPER', 'BVAL_BEST', 'BVAL_LOWER',
               'BVAL_UPPER', 'BVAL_FIX', 'BVAL_FIX_S', 'YCOMP', 'MCOMP',
               'CAT_YMAX', 'PREF_STK', 'PREF_DIP', 'PREF_RKE', 'SHMAX',
               'SHMAX_SIG', 'TRT', 'GMM_TRT', 'DOMAIN', 'CAT_FILE']

fieldType = ['C', 'C', 'C', 'C', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F',
               'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'F', 'C',
               'C', 'F', 'F', 'F', 'F', 'F', 'F', 'C', 'C', 'F', 'C']

fieldSize = [100,  12,  10,  10,   8,   6,   6,   6,   6,   6,   6,   3,   4,
                 4,   4,   4,   4,   8,   8,   8,   6,   6,   6,   6,   6,  70,
                50,   8,   6,   6,   6,   6,   6, 100, 100,   2,  50]

fieldDecimal = [0, 0, 0, 0, 2, 4, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 2, 5, 5, 5, 3, 3, 3,
                3, 3, 0, 0, 3, 2, 2, 2, 2, 2, 0, 0, 0, 0]

# get field data
src_codes = get_field_data(sf, 'CODE', 'str')


#w.field('SRC_NAME','C','100')


outshp = 'test_NSHA18.shp'
'''
build_source_shape(outshp, shapes, src_names, src_codes, zone_class, \
                   rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                   min_rmag, mmax, bval_fix, bval_sig_fix, \
                   ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                   shmax_pref, shmax_sig, trt_new, domains, prefCat)

'''