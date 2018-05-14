import numpy as np
import matplotlib.pyplot as plt
from os import path

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

hxfix = 6. #4.0 # hinge magnitude
def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    #hxfix = 5.5 #4.0 # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hxfix)
    modx_hi = highside(x, hxfix)

    ans1 = modx_lo * (c[0] * x + c[1])
    yarea = c[0] * hxfix + c[1]
    ans2 = modx_hi * (c[2] * (x-hxfix) + yarea)

    return ans1 + ans2

def quadratic_vertex_reg(c, x):
    # set vertex to M8.0
    vertex = 8.0
    xx = x - vertex
        
    return c[0] * xx**2 + c[1]
    
def quadratic_reg(c, x):
    return c[0] * x**2 + c[1]
    
def quadratic_slant_reg(c, x):
    return ((c[0] * x**2) - x) / x # + c[1]
    
def linear_reg(c, x):
    return c[0]*x + c[1]



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

# parse file
nsha_file = path.join('..','data','NSHA18CAT.MS-MW.csv')

cat_nsha = np.genfromtxt(nsha_file,delimiter=',',skip_header=1,dtype=None)

year = []
ms = []
mw = []
msauth = []
idx = []
for x in cat_nsha:
    year.append(int(x[0][0:4]))
    ms.append(x[10])
    mw.append(x[8])
    msauth.append(x[11])
    idx.append(0)

year = np.array(year)
msauth = np.array(msauth)
ms = np.array(ms)
mw = np.array(mw)
idx_src = np.array(idx)

# ignore GA MSmax!
msmidx = np.where((year >= 2014) & (msauth == 'AUST'))
ms = np.delete(ms, msmidx)
mw = np.delete(mw, msmidx)
idx_src = np.delete(idx_src, msmidx)

# get rid of southern ocean event
msmidx = np.where(mw > 7.0)
ms = np.delete(ms, msmidx)
mw = np.delete(mw, msmidx)
idx_src = np.delete(idx_src, msmidx)

ms = np.hstack((ms,[8]))
mw = np.hstack((mw,[8]))


# Mw = np.array([x[3] for x in cat])
# Ml_pre = np.array([x[4] for x in cat])
# Ml_rev = np.array([x[5] for x in cat])
# Mw_ref = np.array([x[6] for x in cat])
#
# idx1 = np.where(Mw_ref=='HG')[0]
# idx2 = np.where(Mw_ref=='TA-SEA')[0]
# idx3 = np.where(Mw_ref=='TA-WA')[0]
# idx4 = np.where(Mw_ref=='other')[0]

############### bilinear
data = odrpack.RealData(ms[ms>=3.5], mw[ms>=3.5])

sx = np.interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(ms[ms>=3.5], mw[ms>=3.5], sx=sx[ms>=3.5], sy=sy[ms>=3.5])

xrng = np.arange(3.5,7.0,step=0.1)

bilin_reg = odrpack.Model(bilinear_reg_free)
odr = odrpack.ODR(data, bilin_reg, beta0=[0.7, 1.0, 1.0, 3.5])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
#out.pprint()

a = out.beta[0]
b = out.beta[1]
c = out.beta[2]
hx = out.beta[3] # x hinge point

yrng = b + a * xrng
yhinge = b + a * hx
idx = xrng > hx
yrng[idx] = c * (xrng[idx]-hx) + yhinge

#hxfix = 5.5
bilin_fix = odrpack.Model(bilinear_reg_fix)
odr = odrpack.ODR(data, bilin_fix, beta0=[0.7, 1., 1.0])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
#print '\nManual\n'
#out.pprint()

af = out.beta[0]
bf = out.beta[1]
cf = out.beta[2]

yrngf = bf + af * xrng
yhinge = bf + af * hxfix
idx = xrng > hxfix
yrngf[idx] = cf * (xrng[idx]-hxfix) + yhinge

############### quadratic vertex

# weight larger MW

quad_reg = odrpack.Model(quadratic_reg)
odr = odrpack.ODR(data, quad_reg, beta0=[1.0, 1.0])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()

print '\n\n'
out.pprint()

aq = out.beta[0]
bq = out.beta[1]
#cq = out.beta[2]
        
quadfit = aq * xrng**2 + bq 
#quadfit= ((aq * xrng**2) - xrng) / xrng + bq
#quadfit = aq*xrng + bq

# write quadratic coefs
coeftxt = 'mw = c1 * ms**2 + c2\n'
coeftxt += str(aq) + '\n'
coeftxt += str(bq) + '\n'
f = open('ms2mw_coefs.txt', 'wb')
f.write(coeftxt)
f.close()

####### linear
#C1 = orthoregress(ms[ms>=3.5], mw[ms>=3.5])
#y1 = C1[0]*xrng + C1[1]
#
#
#
# ################## TA-model
mw_ta = 0.060977193461713*(xrng-7.)**2 + 0.988674360853538*(xrng-7.) + 7.
# ################## Johnston
log_mw_j = 24.66 - 1.083*xrng + 0.192*xrng**2
mw_j = 2*log_mw_j / 3 - 10.7
# ################## Scordilis (2006)
mw_s = np.concatenate((0.67 * xrng[xrng<=6.1]+2.07,
                        0.99 * xrng[xrng>6.1] + 0.08))
# ################## Das etal 2010
mw_das = np.concatenate((0.67 * xrng[xrng<=6.1]+2.12,
                        1.06 * xrng[xrng>6.1] - 0.38))
#mw_das = 1.538*xrng - 2.538
# ################## Youngs (2012)
mw_yon = 2.654 + 0.334*xrng + 0.040*xrng**2

# ################## Di Giacomo (2015)
mw_dig = np.concatenate((0.67 * xrng[xrng<=6.47]+2.13,
                        1.10 * xrng[xrng>6.47] - 0.67))
#
f, ax = plt.subplots(1, 1,figsize=(10,10))

ax.plot([1.5,7.0],[1.5,7.0],'k--',lw=2, label='1:1')
ax.set_xlim([3.5,7.0])
ax.set_ylim([3.5,7.0])
#ax.scatter(ms[idx_src==1],mw[idx_src==1],s= 200, alpha=0.5, marker='o',c='b',label='Ghasemi et al. (2016)')
ax.scatter(ms[idx_src==0],mw[idx_src==0],s= 150, alpha=0.5, marker='o',c='r',label='Data')
# ax.scatter(Ml_pre[idx2],Mw[idx2],s= 70, alpha=0.5, marker=(5,1),c='y',label='TA-SEA')
# ax.scatter(Ml_pre[idx3],Mw[idx3],s= 70, alpha=0.5, marker='^',c='r',label='TA-WA')
# ax.scatter(Ml_pre[idx4],Mw[idx4],s= 70, alpha=0.5, marker='>',c='m',label='Other')
#
# ax.plot(xrng,y1,'g--',lw=3,label='Linear Fit')

ax.plot(xrng,mw_j,'m-',lw=2,label='Johnston (1996)')
ax.plot(xrng,mw_s,'y-',lw=2,label='Scordilis (2006)')
ax.plot(xrng,mw_das,'b-',lw=2,label='Das et al. (2010)')
ax.plot(xrng,mw_yon,'c-',lw=2,label='Youngs (2012)')
ax.plot(xrng,mw_ta,'r-',lw=2,label='Allen (2012; Unpub.)')
ax.plot(xrng,mw_dig,'g-',lw=2,label='Di Giacomo et al (2015)')
#ax.plot(xrng,yrng,'b-',lw=2,label='Automatic-fit')
ax.plot(xrng,quadfit,'k-',lw=2.5,label='Quadratic')
#ax.plot(xrng,yrngf,'b-',lw=2,label='Manual-fit')
#
#
ax.set_xlabel(r'$\mathregular{M_{S}}$', fontsize=22)
ax.set_ylabel(r'$\mathregular{M_{W}}$', fontsize=22)
ax.legend(loc="upper left",ncol=1, scatterpoints=1,fontsize=16)
ax.set_aspect('equal')
ax.grid(which='major')
ax.tick_params(axis='both', labelsize=17)
#
#plt.savefig('MW_MS.png',dpi=300,bbox_inches='tight')
plt.savefig('ms2mw.png',dpi=300,bbox_inches='tight')
plt.show()
plt.close()