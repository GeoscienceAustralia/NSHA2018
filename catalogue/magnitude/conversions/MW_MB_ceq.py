import numpy as np
import matplotlib.pyplot as plt

import scipy.odr.odrpack as odrpack

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x >= hx
    xmod[idx] = 1
    return xmod

def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yhinge = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yhinge)

    return ans1 + ans2

def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    hx = 4.5 #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yarea)

    return ans1 + ans2

from scipy.odr import Model, Data, ODR
from scipy.stats import linregress,norm

def orthoregress(x, y):
    """Perform an Orthogonal Distance Regression on the given data,
    using the same interface as the standard scipy.stats.linregress function.
    Arguments:
    x: x data
    y: y data
    Returns:
    [m, c, nan, nan, nan]
    Uses standard ordinary least squares to estimate the starting parameters
    then uses the scipy.odr interface to the ODRPACK Fortran code to do the
    orthogonal distance calculations.
    """
    linreg = linregress(x, y)
    mod = Model(f)
    dat = Data(x, y)
    od = ODR(dat, mod, beta0=linreg[0:2])
    out = od.run()

    return list(out.beta) + [np.nan, np.nan, np.nan]


def f(p, x):
    """Basic linear regression 'model' for use with ODR"""
    return (p[0] * x) + p[1]


cat_file_TA = '../data/2012.mb_mw_events.csv'
cat_file_KIWI = '../data/mw_kiwi_ml_ga_mb_isc_ms_isc_rm_LP.csv'

cat_ta = np.genfromtxt(cat_file_TA,delimiter=',',skip_header=1,dtype=None)
cat_kiwi = np.genfromtxt(cat_file_KIWI,delimiter=',',skip_header=1,dtype=None)
mb = []
mw = []
idx = []
for x in cat_ta:
    mb.append(x[12])
    mw.append(x[8])
    idx.append(0)
for y in cat_kiwi:
    if not(np.isnan(y[7])):
        mb.append(y[7])
        mw.append(y[3])
        idx.append(1)
mb = np.array(mb)
mw = np.array(mw)
idx_src = np.array(idx)


# Mw = np.array([x[3] for x in cat])
# Ml_pre = np.array([x[4] for x in cat])
# Ml_rev = np.array([x[5] for x in cat])
# Mw_ref = np.array([x[6] for x in cat])
#
# idx1 = np.where(Mw_ref=='HG')[0]
# idx2 = np.where(Mw_ref=='TA-SEA')[0]
# idx3 = np.where(Mw_ref=='TA-WA')[0]
# idx4 = np.where(Mw_ref=='other')[0]
############### bilinear auto
data = odrpack.RealData(mb[mb>=3.5], mw[mb>=3.5])
xrng = np.arange(3.5,6.0,step=0.01)

bilin_reg = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data, bilin_reg, beta0=[0.7, 1.0, 1.0, 3.5])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
out.pprint()

a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3] # x hinge point

yrng = b + a * xrng
yhinge = b + a * hx
idx = xrng > hx
yrng[idx] = c * (xrng[idx]-hx) + yhinge

##################################bilinear fix
hxfix = 4.75
bilin_fix = odrpack.Model(bilinear_reg_fix)
odr = odrpack.ODR(data, bilin_fix, beta0=[0.7, 1., 1.0])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
out.pprint()

af = out.beta[0]
bf = out.beta[1]
cf = out.beta[2]

yrngf = bf + af * xrng
yhinge = bf + af * hxfix
idx = xrng > hxfix
yrngf[idx] = cf * (xrng[idx]-hxfix) + yhinge

####### linear
C1 = orthoregress(mb[mb>=3.5], mw[mb>=3.5])
y1 = C1[0]*xrng + C1[1]
#
#
#
# ################## TA-model
mw_ta = np.concatenate((0.686170473456368 * xrng[xrng<=5.2]+1.152064114214877,
                        0.686170473456368 * xrng[xrng>5.2] + 0.982366026214298 * (xrng[xrng>5.2] - 5.2) + 1.152064114214877))
# ################## Johnson1996
log_mw_J = 18.28 + 0.679*xrng + 0.077*xrng**2
mw_j = (2*log_mw_J/3.) - 10.7

# ################## Scordilis (2006)
mw_s = 0.85*xrng + 1.03

# ################## Sonley & Atkinson 2005
mw_sa = 1.03*xrng - 0.61

# ################## Das etal 2010
mw_das = 0.65*xrng + 1.65

# ################## Youngs (2012)
mw_yon = xrng - 0.42

# ################## Di Giacomo (2015)
mw_dig = 1.38*xrng - 1.79


#
f, ax = plt.subplots(1, 1,figsize=(10,10))

ax.plot([1.5,6.0],[1.5,6.0],'--', c='0.4',lw=2.0, label='1:1')
ax.set_xlim([3.0,6.0])
ax.set_ylim([3.0,6.0])
ax.scatter(mb[idx_src==1],mw[idx_src==1],s= 200, alpha=0.5, marker='o',c='b',label='Ghasemi et al.')
ax.scatter(mb[idx_src==0],mw[idx_src==0],s= 200, alpha=0.5, marker='^',c='r',label='Other')
# ax.scatter(Ml_pre[idx2],Mw[idx2],s= 70, alpha=0.5, marker=(5,1),c='y',label='TA-SEA')
# ax.scatter(Ml_pre[idx3],Mw[idx3],s= 70, alpha=0.5, marker='^',c='r',label='TA-WA')
# ax.scatter(Ml_pre[idx4],Mw[idx4],s= 70, alpha=0.5, marker='>',c='m',label='Other')
#
ax.plot(xrng,mw_j,'m-',lw=2,label='Johnston (1996)')
ax.plot(xrng,mw_sa,'b-',lw=2,label='Sonley & Atkinson (2005)')
ax.plot(xrng,mw_s,'y-',lw=2,label='Scordilis (2006)')
ax.plot(xrng,mw_das,'b--',lw=2,label='Das et al. (2010)')
ax.plot(xrng,mw_yon,'c-',lw=2,label='Youngs (2012)')
ax.plot(xrng,mw_ta,'r-',lw=2,label='Allen (2012)')
ax.plot(xrng,mw_dig,'g-',lw=2,label='Di Giacomo et al. (2015)')
# ax.plot(xrng,yrng,'b-',lw=2,label='Automatic-fit')
ax.plot(xrng,yrngf,'k-',lw=2,label='Fixed Hinge')
ax.plot(xrng,y1,'k--',lw=3,label='Linear')

#
#
#
ax.set_xlabel('mb', fontsize=25)
ax.set_ylabel('MW', fontsize=25)
ax.legend(loc="upper left",ncol=1,scatterpoints=1,fontsize=15)
ax.set_aspect('equal')
ax.grid(which='major')
ax.tick_params(axis='both', labelsize=21)
#
plt.savefig('MW_MB.TA.png',dpi=300,bbox_inches='tight')
plt.close()