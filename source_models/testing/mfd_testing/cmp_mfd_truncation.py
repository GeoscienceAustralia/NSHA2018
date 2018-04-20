# gets incremental earthquake rates for openquake source files
def get_oq_incrementalMFD(beta, N0, mmin, mmax, binwid):
    from numpy import arange, exp

    mrange = arange(mmin, mmax, binwid)

    betacurve = N0 * exp(-beta  *mrange) * (1 - exp(-beta * (mmax - mrange))) \
                / (1 - exp(-beta * mmax))

    return betacurve, mrange

def bval2beta(bval):
    from numpy import log
    return log(10**bval)


from openquake.hazardlib.mfd.truncated_gr import TruncatedGRMFD
#from oq_tools import bval2beta, get_oq_incrementalMFD
import matplotlib.pyplot as plt
from numpy import array, log10

min_mag = 4.5 # magnitude of earthquake
max_mag = 7.5
a_val = log10(4e4) # random number
N0 = 10**a_val
b_val = 1.2
beta = bval2beta(b_val)
bin_width = [0.1] #, 0.01]

fig = plt.figure(1, figsize=(12,6))

# loop through bin width
for i, bw in enumerate(bin_width):
    plt.subplot(1,2,i+1)
    
    # get oq rates
    gr = TruncatedGRMFD(min_mag, max_mag, bw, a_val, b_val)
    
    oq_rates = gr.get_annual_occurrence_rates()
    # get cum OQ ratescum_oq_rates = []
    
    
    # get cummulative frisk rates (for OQ)
    cum_frisk_rates, frisk_mags = get_oq_incrementalMFD(beta, N0, min_mag, max_mag, bw)
    
    # get rates
    frisk_rates = []
    for i in range(0, len(cum_frisk_rates)):
        if i < len(cum_frisk_rates)-1:
            frisk_rates.append(cum_frisk_rates[i] - cum_frisk_rates[i+1])
        else:
            frisk_rates.append(cum_frisk_rates[i])
    
    
    # plt rates
    plt.semilogy(array(oq_rates)[:,0], array(oq_rates)[:,1], 'r-', lw=2.0)
    #plt.semilogy(frisk_mags+0.5*bin_width, frisk_rates, '-', c='lightblue', lw=2.0)
    plt.semilogy(frisk_mags+0.5*bw, frisk_rates, '-', c='b', lw=2.0)

plt.show()

#sum(array(gr.get_annual_occurrence_rates())[:,1])