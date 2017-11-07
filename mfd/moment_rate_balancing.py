"""Code to test moment balancing of fault sources.
This takes source models and compared the resulting mfds, before
testing that the moment balance is met. 
"""

import os, sys
import numpy
import hmtk
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
#from openquake.commonlib.nrml import read

def generate_mfd(source_model_file):
    parser = nrmlSourceModelParser(source_model_file)
    source_model_name = source_model_file.split('/')[-1].rstrip('.xml')
    source_model = parser.read_file(source_model_name)
#    source_model = read(source_model_file)
    print source_model_name
    print source_model
    for source in source_model:
        print source, type(source)
        print source.mfd
        occurrence_rates = source.mfd.get_annual_occurrence_rates()
        minmag = source.mfd._get_min_mag_and_num_bins()
        print minmag
        if str(source.mfd) == "<TruncatedGRMFD>":
            #Theoretical moment rate
            moment_rate_pure = source.mfd._get_total_moment_rate()
        elif str(source.mfd) == "<YoungsCoppersmith1985MFD>":
            # Have to calculate the total moment rate ourselves
            moment_rate_pure = None
        # Now calculate actual discretised moment rate
        moment_rate = 0
        for pair in occurrence_rates:
            print pair
            moment = numpy.power(10, (1.5*(pair[0])+16.05))
            moment = moment/1e7 # Nm
            inc_moment_rate = moment*pair[1]
            moment_rate += inc_moment_rate
        print moment_rate_pure, moment_rate

if __name__ == "__main__":
    path = '../shapefile2nrml/'
    GR_nrml_source_model = \
        os.path.join(path, 'Adelaide_test/Adelaide_faults.xml')
    CE_nrml_source_model = \
        os.path.join(path, 'Adelaide_test_CE/Adelaide_faults.xml')

    GR_mfds = generate_mfd(GR_nrml_source_model)
    CE_mfds = generate_mfd(CE_nrml_source_model)
