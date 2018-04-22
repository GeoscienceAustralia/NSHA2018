import numpy as np
import matplotlib.pyplot as plt

import scipy.odr.odrpack as odrpack
from collections import OrderedDict
import pdb
from misc_tools import checkfloat
from os import path

########
#hxfix = 4.
########

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

hxfix = 4.5 #4.0 # hinge magnitude
def bilinear_reg_fix(c, x):
    from numpy import zeros_like
    hxfix = 4.5 #4.0 # hinge magnitude
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

def f(B, x):
    return B[0]*x**2+B[1]*x+B[2]

####################################################################
# parse file
nsha_file = path.join('..','..','data','NSHA18CAT.ML-MW.csv')

lines = open(nsha_file).readlines()[1:]

ml = []
mlr = []
mw = []
mwref = []
idx = []

for line in lines:
    x = line.strip().split(',')
    if np.isnan(checkfloat(x[17])):
        ml.append(float(x[14]))
        print checkfloat(x[17])
    else:
        ml.append(float(x[17]))
    mw.append(float(x[8]))
    mwref.append(x[9])
    idx.append(0)

ml = np.array(ml)
mw = np.array(mw)
mwref = np.array(mwref)
idx_src = np.array(idx)

# delete events
didx = np.where((ml < 3.75) & (mwref=='Ghasemi et al (2016)'))[0]
ml = np.delete(ml, didx)
mw = np.delete(mw, didx)
mwref = np.delete(mwref, didx)

didx = np.where(mw > 6.25)[0]
print ml[didx]
ml = np.delete(ml, didx)
mw = np.delete(mw, didx)
mwref = np.delete(mwref, didx)



'''mw = np.array([x[3] for x in cat])
ml = np.array([x[4] for x in cat])
ml = np.array([x[5] for x in cat])
mw_ref = np.array([x[6] for x in cat])
'''
# # not include mw<4.0 for HG estimates
# idx_exclude = np.where((mw_ref=='HG') & (mw<4.0))[0]
# mw = mw[~idx_exclude]
# ml = ml[~idx_exclude]
# ml = ml[~idx_exclude]
# mw_ref = mw_ref[~idx_exclude]


idx1 = np.where(mwref=='Ghasemi et al (2016)')[0]
idx2 = np.where(mwref=='Allen (2012)')[0]
idx3 = np.where(mwref=='Allen et al. (2006)')[0]
idx4 = []
for i, mr in enumerate(mwref):
    if mr != 'Ghasemi et al (2016)' and mr != 'Allen (2012)' and mr != 'Allen et al. (2006)':
        idx4.append(i)

############### bilinear auto
#data = odrpack.RealData(ml, mw)
data = odrpack.RealData(ml, mw)

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
dat = Data(ml, mw)
co = np.polynomial.polynomial.polyfit(ml, mw,2)
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
# pw, cov = curve_fit(two_lines, ml, mw, pw0)
# yrngf_mix = two_lines(xrng,*pw)
################## Swiss
yrng_swiss = np.concatenate((0.594 * xrng[xrng<=2]+0.985,
                             1.327+0.253*xrng[(xrng>2)&(xrng<=4)]+0.085*xrng[(xrng>2)&(xrng<=4)]**2,
                             xrng[xrng>4]-0.3))
# yrng_swiss = 1.02 + 0.472 * xrng + 0.0491 * xrng ** 2

##### Ross et al. 2016

yrng_ross = 0.754*xrng[xrng<=4.0] + 0.88


f, ax = plt.subplots(1, 1, figsize=(10,10))
ax.plot([2.5,7.0],[2.5,7.0],'k--',lw=1,label='1:1')
ax.set_xlim([2.5,7.0])
ax.set_ylim([2.5,7.0])
ax.scatter(ml[idx3],mw[idx3],s= 70, alpha=0.5, marker='^',c='r',label='Allen et al. (2006)')
ax.scatter(ml[idx2],mw[idx2],s= 70, alpha=0.5, marker=(5,1),c='y',label='Allen (2012)')
ax.scatter(ml[idx1],mw[idx1],s= 70, alpha=0.5, marker='o',c='b',label='Ghasemi et al. (2017)')
ax.scatter(ml[idx4],mw[idx4],s= 70, alpha=0.5, marker='>',c='m',label='Other')

ax.plot(xrng,yrng_swiss,'r-',lw=2,label='Goertz-Allmann et al. (2011)')
ax.plot(xrng[xrng<=4.0],yrng_ross,'b-',lw=2,label='Ross et al. (2016)')
ax.plot(xrng,yrng,'b-',lw=2,label='Automatic-fit')
ax.plot(xrng,yrngf,'g-',lw=2,label='Bi-linear')
ax.plot(xrng,yrng_poly,'k-',lw=2,label='Polynomial')

# ax.plot(xrng,yrngf_mix,'m-',lw=4)

handles, labels = ax.get_legend_handles_labels()

labels = [labels[0],labels[7],labels[6],labels[5],labels[8],labels[3],labels[4],labels[1],labels[2]]
handles = [handles[0],handles[7],handles[6],handles[5],handles[8],handles[3],handles[4],handles[1],handles[2]]

#leg = ax.legend(handles,labels,loc="upper left",ncol=1,scatterpoints=1)

# make pretty
ax.set_xlabel(r'$\mathregular{M_{L}}$', fontsize=22)
ax.set_ylabel(r'$\mathregular{M_{W}}$', fontsize=22)
ax.legend(loc="upper left",ncol=1, scatterpoints=1,fontsize=16)
ax.set_aspect('equal')
ax.grid(which='major')
ax.tick_params(axis='both', labelsize=17)
#
plt.savefig('ml2mw.png',dpi=300,bbox_inches='tight')
plt.show()
plt.close()