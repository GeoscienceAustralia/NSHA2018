# -*- coding: utf-8 -*-
"""
Created on Wed Feb 07 16:11:14 2018

@author: u56903
"""

from parsers import parse_iscgem
from writers import iscgem2hmtk_csv
from os import path

# parse catalogue
iscgemDict = parse_iscgem('data/isc-gem-cat.csv')

# convert to hmtk
iscgem2hmtk_csv(iscgemDict, path.join('data', 'ISC-GEM_V4_hmtk_full.csv'))
 