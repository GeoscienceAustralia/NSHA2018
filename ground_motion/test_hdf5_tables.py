from calc_oq_gmpes import hdf5_gsim, interface_gsims
from numpy import logspace, sqrt, array, exp
import matplotlib.pyplot as plt
from os import path, getcwd
from scipy.constants import g
import matplotlib as mpl
mpl.style.use('classic')


mag  = 7.0
dep = 50.
ztor = 40. # guess
rake = -90. # USGS CMT
dip  = 60.

# set site details
vs30 = 760.
#rjb = logspace(1,2.7,10)
rjb = array([100, 800])

# at this distance, assume all distances are equivalent
rrup = rjb #sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
rhypo = rjb

fig = plt.figure(1, figsize=(18, 9))
wdir = getcwd()

for i, r in enumerate(rjb):

    # get ground motion estimates from GMPEs
    #for i in range(0, len(rrup)):
    #Yea97imt, AB03imt, AB03CISimt, Gea05imt, Zea06imt, Zea06CISimt, MP10imt, AA13imt, Aea15imt, Zea16imt \
    Yea97imt, AB03imt, Zea06imt, Zea06CISimt, AM09imt, MP10imt, GA14imt, GA14CISimt, Aea15imt \
        = interface_gsims(mag, dep, ztor, dip, rake, r, r, vs30)
    
    # plot Garcia
    plt.subplot(1,2,i+1)    
    hdf5file = path.join(wdir, 'hdf5_tables', 'MegawatiPan2010.vs760.h30.hdf5')
    MPhdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, r, r, r, vs30, hdf5file)
    # M6
    
    plt.loglog(array(MPhdf5imt['per']), exp(MPhdf5imt['sa']), 'r-', lw=2, label='GSIM Table')
    plt.loglog(array(MP10imt['per']), exp(MP10imt['sa']), 'b-', lw=2, label='Parametric')
    plt.loglog(array(Aea15imt['per']), exp(Aea15imt['sa']), 'k--', lw=2, label='Extrap Host (AB15)')
    
    # plot from table
    '''
    if i == 1:
       tabPer = Geahdf5imt['per'][::-1]
       tabSA = array([2.477, 2.541, 2.633, 2.667, 2.503, 2.348, 2.159, 2.019, 1.748, 1.560, 1.262, 1.030, 0.697, 0.445, 0.250, -0.294])
       plt.loglog(tabPer, (10**tabSA)/(g*100.), 'g--', lw=2, label='ASCII Table')
    '''
    plt.legend()
    
    plt.xlabel('Period (sec)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.title(' '.join(('Megawatti & Pan (2005); MW =', str(mag) + '; Rrup =',str('%0.1f' % r), 'km; h =',str('%0.0f' % dep), 'km')))
    plt.grid(which='both')
    plt.xlim([0.01, 12.])
    
    """
    plt.subplot(122)    
    hdf5file = path.join(wdir, 'hdf5', 'A.15_SP15_adjusted_spec.hdf5')
    SPhdf5imt = hdf5_gsim(mag, dep, ztor, dip, rake, rrup[1], rjb[1], vs30, hdf5file)
    tabtxt = '-0.278 0.011 0.360 0.541 0.773 1.086 1.298 1.584 1.776 2.027 2.154 2.304 2.391 2.488 2.595 2.707 2.754 2.769 2.752 2.710 2.671 2.624 2.481'
    # M5
    #tabtxt = '-1.374 -0.966 -0.596 -0.418 -0.177 0.190 0.463 0.850 1.113 1.455 1.625 1.824 1.938 2.066 2.211 2.371 2.450 2.502 2.500 2.469 2.431 2.385 2.230'
    tabvals = (10**array([float(x) for x in tabtxt.split()]) / (100 * g))[::-1]
    
    plt.loglog(array(SP16imt['per']), exp(SP16imt['sa']), 'b-', lw=2, label='Parametric')
    plt.loglog(array(SPhdf5imt['per']), exp(SPhdf5imt['sa']), 'r-', lw=2, label='GSIM Table')
    plt.loglog(array(SPhdf5imt['per']), tabvals, 'g--', lw=2, label='Tables')
    plt.legend()
    
    plt.xlabel('Period (sec)')
    plt.ylabel('Spectral Acceleration (g)')
    plt.title(' '.join(('Shahjouei & Pezeshk (2016); MW =', str(mag) + '; Rrup =',str('%0.1f' % rrup[1]), 'km')))
    plt.grid(which='both')
    plt.xlim([0.01, 10.])
    """
plt.savefig('hdf5_test_M'+str(mag)+'_vs'+str('%0.0f' % vs30)+'.png', fmt='png', bbox_inches='tight')

plt.show()

