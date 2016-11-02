"""Example code showing differences in input vs output 
moment rate for YC1985 distirbution
"""

from openquake.hazardlib.mfd import YoungsCoppersmith1985MFD
import numpy as np
from matplotlib import pyplot

def YC_moment_rates(characteristic_mag, b_value, min_mag, moment_rate, bin_width):
    """
    Generate Youngs and Coppersmith MFD using OpenQuake, then convert
    to incrementalMFD type. This allows us to balance the seismic moment to 
    be exactly equal to the total moment rate. We do this assuming the GR
    part of the YC distribution starts at magnitude 0.01.
    """

    mfd = YoungsCoppersmith1985MFD.from_total_moment_rate(min_mag=0.01,
                                                          b_val=float(b_value),
                                                          char_mag=characteristic_mag, 
                                                          total_moment_rate=moment_rate,
                                                          bin_width=float(bin_width))
    mags,rates=zip(*mfd.get_annual_occurrence_rates())
    #print mags, rates
    # calcualate total moment rate and rescale rates if
    # necessary to meet total input rate
    total_moment_rate = 0
    for i in range(len(mags)):
        moment = np.power(10, (1.5*mags[i]+16.05))
        moment = moment/1e7 #Nm
        inc_moment_rate = moment*rates[i]
        total_moment_rate += inc_moment_rate
    moment_error = (total_moment_rate - moment_rate)/moment_rate
    initial_moment_error = moment_error
    print 'Input moment rate ', moment_rate
    print 'MFD moment rate', total_moment_rate
    print 'Relative moment rate error', moment_error

    # Rescale rates
#    print rates
    rates = rates/(1+moment_error)
#    print rates
    # check rates sum as expected
    total_moment_rate = 0
    for i in range(len(mags)):
        moment = np.power(10, (1.5*mags[i]+16.05))
        moment = moment/1e7 #Nm
        inc_moment_rate = moment*rates[i]
        total_moment_rate += inc_moment_rate
    moment_error = (total_moment_rate - moment_rate)/moment_rate
    print 'Final moment rate error after rescaling',  moment_error

    # Now trim the distribution to just above min_mag
#    print mags, rates
    mags = np.array(mags)
    rates = rates[np.where(mags >= float(min_mag))]
    mags = mags[np.where(mags >= float(min_mag))]
    return initial_moment_error


if __name__ == "__main__":
    """
    """
    char_mag = 7.2 # Smallest errors on intervals of 0.05, e.g. 7.25, largest errors for 7.2, 7.3
    b_value = 1.0
    min_mag = 4.5 # This has no effect as in code above we calculate GR part of MFD down to 0.01
    moment_rate = 1.3e15
    bin_width = 0.1 # Reducing bin width to 0.001 reduces error
    mag_list = []
    error_list = []
    for mag in np.arange(7.0, 7.4, 0.01):
        char_mag = mag
        initial_moment_error = YC_moment_rates(char_mag, b_value, min_mag,
                                               moment_rate, bin_width)
        mag_list.append(mag)
        error_list.append(initial_moment_error)
    pyplot.scatter(mag_list, error_list)
    pyplot.xlabel('Characteristic magnitude')
    pyplot.ylabel('(MFD moment rate - initial moment rate)/initial moment rate')
    pyplot.savefig('moment_rate_relative_error.png')
