"""Code to test moment balancing of fault sources.
This takes source models and compared the resulting mfds, before
testing that the moment balance is met. 
"""

import os, sys
import numpy
import hmtk
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser

def generate_mfd(source_model_file):
    parser = nrmlSourceModelParser(source_model_file)
    source_model_name = source_model_file.split('/')[-1].rstrip('.xml')
    source_model = parser.read_file(source_model_name)
    print source_model_name
    for source in source_model:
        print source
        print source.mfd
        occurrence_rates = source.mfd.get_annual_occurrence_rates()
        print occurrence_rates
        minmag = source.mfd._get_min_mag_and_num_bins()
        print minmag
        moment_rate = source.mfd._get_total_moment_rate()
        print moment_rate

if __name__ == "__main__":
    path = '../shapefile2nrml/'
    GR_nrml_source_model = \
        os.path.join(path, 'Adelaide_test/Adelaide_faults.xml')
    CE_nrml_source_model = \
        os.path.join(path, 'Adelaide_test_CE/Adelaide_faults.xml')

    GR_mfds = generate_mfd(GR_nrml_source_model)
    CE_mfds = generate_mfd(CE_nrml_source_model)
