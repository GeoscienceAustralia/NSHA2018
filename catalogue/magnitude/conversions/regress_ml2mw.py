import numpy as np
import matplotlib.pyplot as plt

import scipy.odr.odrpack as odrpack
from collections import OrderedDict
import pdb

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

def f(B, x):
    return B[0]*x**2+B[1]*x+B[2]

# cat_file = '../data/MW_ML_all_test_2.csv'
cat_file = '../data/MW_ML_all_test_3_fill_Ml_rev.csv'

cat = np.genfromtxt(cat_file,delimiter=',',skip_header=1,dtype=None)

Mw = np.array([x[3] for x in cat])
Ml_pre = np.array([x[4] for x in cat])
Ml_rev = np.array([x[5] for x in cat])
Mw_ref = np.array([x[6] for x in cat])

# # not include Mw<4.0 for HG estimates
# idx_exclude = np.where((Mw_ref=='HG') & (Mw<4.0))[0]
# Mw = Mw[~idx_exclude]
# Ml_pre = Ml_pre[~idx_exclude]
# Ml_rev = Ml_rev[~idx_exclude]
# Mw_ref = Mw_ref[~idx_exclude]


idx1 = np.where(Mw_ref=='HG')[0]
idx2 = np.where(Mw_ref=='TA-SEA')[0]
idx3 = np.where(Mw_ref=='TA-WA')[0]
idx4 = np.where(Mw_ref=='other')[0]

############### bilinear auto
#data = odrpack.RealData(Ml_pre, Mw)
data = odrpack.RealData(Ml_rev, Mw)

xrng = np.arange(1.5,6.0,step=0.1)

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

###############bilinear fix
hxfix = 4.5
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
# pdb.set_trace()
############## polynomial
from scipy.odr import Model, Data, ODR
from scipy.stats import linregress,norm
mod = Model(f)
dat = Data(Ml_rev, Mw)
co = np.polynomial.polynomial.polyfit(Ml_rev, Mw,2)
od = ODR(dat, mod, beta0=[co[2],co[1],co[0]])
out = od.run()

yrng_poly = out.beta[0]*xrng**2+out.beta[1]*xrng+out.beta[2]

out_poly = out

# ############## polynomial + linear
# def mix_reg_fix(c, x):
#     from numpy import zeros_like
#     hx = 4.5 #4.0 # hinge magnitude
#     ans2 = zeros_like(x)
#     ans1 = zeros_like(x)
#
#     modx_lo = lowside(x, hx)
#     modx_hi = highside(x, hx)
#
#     ans1 = modx_lo * (c[0] * x ** 2 + c[1] * x + c[2])
#     yarea = c[0] * hx ** 2 + c[1] * hx +c[2]
#     ans2 = modx_hi * (c[3] * (x-hx) + yarea)
#
#     return ans1 + ans2
#
# polylin_fix = odrpack.Model(mix_reg_fix)
# odr = odrpack.ODR(data, polylin_fix, beta0=[0.1, 1., 1.0, 1.0])
#
#
# odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
# out = odr.run()
# out.pprint()
# xrng = np.arange(1.5,6.0,step=0.1)
# af = out.beta[0]
# bf = out.beta[1]
# cf = out.beta[2]
# df = out.beta[3]
#
# yrngf_mix =  af * xrng ** 2 + bf * xrng + cf
# yhinge_mix = af * hxfix ** 2 + bf * hxfix + cf
# idx = xrng > hxfix
# yrngf_mix[idx] = df * (xrng[idx]-hxfix) + yhinge_mix

# from scipy.optimize import curve_fit
#
# def two_lines(x, a, b, c, d, e):
#     one = a*x**2 + b*x + c
#     two = d*x + e
#     return np.maximum(one, two)
# pw0 = (1, 1, 1, 1,1) # a guess for slope, intercept, slope, intercept
# pw, cov = curve_fit(two_lines, Ml_rev, Mw, pw0)
# yrngf_mix = two_lines(xrng,*pw)
################## Swiss
yrng_swiss = np.concatenate((0.594 * xrng[xrng<=2]+0.985,
                             1.327+0.253*xrng[(xrng>2)&(xrng<=4)]+0.085*xrng[(xrng>2)&(xrng<=4)]**2,
                             xrng[xrng>4]-0.3))
# yrng_swiss = 1.02 + 0.472 * xrng + 0.0491 * xrng ** 2

##### Ross et al. 2016

yrng_ross = 0.754*xrng[xrng<=4.0] + 0.88


f, ax = plt.subplots(1, 1)
f.set_size_inches(10, 10)
ax.plot([1.5,6.0],[1.5,6.0],'k--',lw=1,label='1:1')
ax.set_xlim([1.5,6.0])
ax.set_ylim([1.5,6.0])
ax.scatter(Ml_pre[idx1],Mw[idx1],s= 70, alpha=0.5, marker='o',c='b',label='Ghasemi et al. (2017)')
ax.scatter(Ml_pre[idx2],Mw[idx2],s= 70, alpha=0.5, marker=(5,1),c='y',label='Allen (2012)')
ax.scatter(Ml_pre[idx3],Mw[idx3],s= 70, alpha=0.5, marker='^',c='r',label='Allen et al. (2006)')
ax.scatter(Ml_pre[idx4],Mw[idx4],s= 70, alpha=0.5, marker='>',c='m',label='Other')

ax.plot(xrng,yrng_swiss,'r-',lw=2,label='Goertz-Allmann et al. (2011)')
ax.plot(xrng[xrng<=4.0],yrng_ross,'b-',lw=2,label='Ross et al. (2016)')
#ax.plot(xrng,yrng,'b-',lw=2,label='Automatic-fit')
ax.plot(xrng,yrngf,'g-',lw=2,label='Bi-linear')
ax.plot(xrng,yrng_poly,'k-',lw=2,label='Polynomial')

# ax.plot(xrng,yrngf_mix,'m-',lw=4)

ax.set_xlabel('ML (revised)',fontsize=16)
ax.set_ylabel('MW',fontsize=16)
ax.grid(which='both')
handles, labels = ax.get_legend_handles_labels()
# pdb.set_trace()
labels = [labels[0],labels[7],labels[6],labels[5],labels[8],labels[3],labels[4],labels[1],labels[2]]
handles = [handles[0],handles[7],handles[6],handles[5],handles[8],handles[3],handles[4],handles[1],handles[2]]
# by_label = OrderedDict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys())
# sort both labels and handles by labels
#labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
leg = ax.legend(handles,labels,loc="upper left",ncol=1,scatterpoints=1)
#leg.get_frame().set_alpha(0.5)
#ax.legend(loc="upper left",ncol=1,numpoints=1)
ax.set_aspect('equal')
plt.show()
# plt.savefig('test.png',bbox_inches='tight',dpi=600)
# plt.close()