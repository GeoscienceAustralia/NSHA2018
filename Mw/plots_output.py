import pdb
import numpy as np
import matplotlib.pyplot as plt
from operator import sub
import matplotlib.patches as mpatch
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap

cat_file = './catalogue.txt'
kiwi_result_dir = '../RESULTS/'
cat = np.loadtxt(cat_file)
id_eqs = []
Mw_ref = []
Ml_ga = []
Ms_ga = []
Mw_sim = []
M0_sim = []
M0_err_sim = []
misfit_sim = []

##m = Basemap(projection='cyl',
##            resolution=None,llcrnrlon=110.,llcrnrlat=-45.,urcrnrlon=155.,urcrnrlat= -10.)
##m.bluemarble()
##plt.show()


def parse_kiwi_output(kiwi_output_file,kiwi_bootstrap_file):
    with open(kiwi_output_file) as f:
        lines = f.readlines()
    Mw = float(lines[23].split(' ')[-1])
    mis_fit = float(lines[29].split(' ')[-1])
    no_traces = int(lines[39].split(' ')[6])
    with open(kiwi_bootstrap_file) as f2:
        lines2 = f2.readlines()
    str1 = lines2[2].split(' ')
    M0 = float(str1[2])
    M0_lower = str1[3].split(',')[0]
    M0_lower = float(M0_lower.split('[')[1])
    M0_upper = float(str1[4].split(']')[0])
    M0_error = [M0_lower,M0_upper]
    return Mw,mis_fit,no_traces,M0,M0_error


for e in cat:
    eve_id = 'auto_eve_' + str(int(e[0]))
    # parse the kiwi outpit file
    kiwi_output_file = kiwi_result_dir + eve_id + '/step1.earthquakeinfo.dat'
    kiwi_bootstrap_file = kiwi_result_dir + eve_id + '/bootstrap.dat'
    Mw_eve,misfit_eve,no_tr_eve,M0_eve,M0_error_eve = parse_kiwi_output(kiwi_output_file,\
                                                    kiwi_bootstrap_file)
    # generate lists for simulation/reference values
    id_eqs.append(eve_id)
    Mw_ref.append(e[11])
    Ml_ga.append(e[12])
    Ms_ga.append(e[13])
    Mw_sim.append(Mw_eve)
    M0_sim.append(M0_eve)
    M0_err_sim.append(M0_error_eve)
    misfit_sim.append(misfit_eve)
# PLOTS
# plot 1: Res: Mw_ref-Mw_sim
res1 = map(sub,Mw_ref,Mw_sim)
res2 = map(sub,Mw_ref,Ml_ga)
res3 = map(sub,Mw_ref,Ms_ga)
res1 = [round(x,2) for x in res1]
res2 = [round(x,2) for x in res2]
res3 = [round(x,2) for x in res3]
fig = plt.figure(figsize=(10,10))
ax = plt.axes()
xlabels = id_eqs
ax.scatter(range(len(res1)),res1,80,marker='o',color='k',alpha=0.5)
ax.scatter(range(len(res1)),res2,80,marker='^',color='r',alpha=0.5)
ax.scatter(range(len(res1)),res3,80,marker='v',color='b',alpha=0.5)
ax.set_xticks(range(len(res1)))
ax.set_xticklabels(xlabels,fontsize=12)
ax.set_xlim([-0.5,len(res1)-0.5])
ax.plot([-0.5,len(res1)-0.5],[0,0],'k-')
ax.set_ylim([-1.1,1.1])
ax.set_yticks(np.arange(-1.1,1.2,0.1))
ax.grid()
ax.legend(['Baseline','This study','GA-Ml','GA-Ms'])
ax.set_ylabel('Error : Mw(ref)-Mw(sim)') 

# plot 2: bar plot
x = [0.1,0.2,0.3,0.4]
y1 = len(np.where(np.abs(res1)<=0.1)[0])
y2 = len(np.where((np.abs(res1)<=0.2)&(np.abs(res1)>0.1))[0])
y3 = len(np.where((np.abs(res1)<=0.3)&(np.abs(res1)>0.2))[0])
y4 = len(np.where((np.abs(res1)>0.3))[0])
y_this = np.array([y1,y2,y3,y4])

y1 = len(np.where(np.abs(res1[6:])<=0.1)[0])
y2 = len(np.where((np.abs(res1[6:])<=0.2)&(np.abs(res1[6:])>0.1))[0])
y3 = len(np.where((np.abs(res1[6:])<=0.3)&(np.abs(res1[6:])>0.2))[0])
y4 = len(np.where((np.abs(res1[6:])>0.3))[0])
y_this_available_ms = np.array([y1,y2,y3,y4])

y1 = len(np.where(np.abs(res2)<=0.1)[0])
y2 = len(np.where((np.abs(res2)<=0.2)&(np.abs(res2)>0.1))[0])
y3 = len(np.where((np.abs(res2)<=0.3)&(np.abs(res2)>0.2))[0])
y4 = len(np.where((np.abs(res2)>0.3))[0])
y_ml_ga = np.array([y1,y2,y3,y4])

y1 = len(np.where(np.abs(res3)<=0.1)[0])
y2 = len(np.where((np.abs(res3)<=0.2)&(np.abs(res3)>0.1))[0])
y3 = len(np.where((np.abs(res3)<=0.3)&(np.abs(res3)>0.2))[0])
y4 = len(np.where((np.abs(res3)>0.3))[0])
y_ms_ga = np.array([y1,y2,y3,y4])

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax1.bar(x,y_this,facecolor='#9999ff',edgecolor='white',width=0.05,alpha=0.5)
for x1,y1 in zip(x,y_this):
    ax1.text(x1+0.05,y1+0.05,'%d'%y1)

ax1.bar(x,-y_ml_ga,facecolor='#ff9999',edgecolor='white',width=0.05,alpha=0.5)
for x1,y1 in zip(x,y_ml_ga):
    ax1.text(x1+0.05,-y1-0.05,'%d'%y1)

ax1.set_xticks(np.array(x)+0.05)
ax1.set_xticklabels(('<0.1','0.1-0.2','0.2-0.3','>0.3'),fontsize=12)
ax1.set_xlim([0.05,0.5])
ax1.set_ylim([np.min([-y_ml_ga,-y_ms_ga])-1,np.max(y_this)+1])
ax1.set_xlabel('Error range')
ax1.set_ylabel('Frequency')
ax1.set_title('Mw (this study) - Ml (GA)')

ax2 = fig.add_subplot(1,2,2)
ax2.bar(x,y_this_available_ms,facecolor='#9999ff',edgecolor='white',width=0.05,alpha=0.5)
for x1,y1 in zip(x,y_this_available_ms):
    ax2.text(x1+0.05,y1+0.05,'%d'%y1)
ax2.bar(x,-y_ms_ga,facecolor='#ff9999',edgecolor='white',width=0.05,alpha=0.5)
for x1,y1 in zip(x,y_ms_ga):
    ax2.text(x1+0.05,-y1-0.05,'%d'%y1)

ax2.set_xticks(np.array(x)+0.05)
ax2.set_xticklabels(('<0.1','0.1-0.2','0.2-0.3','>0.3'),fontsize=12)
ax2.set_xlim([0.05,0.5])
ax2.set_ylim([np.min([-y_ml_ga,-y_ms_ga])-1,np.max(y_this)+1])
ax2.set_xlabel('Error range')
ax2.set_title('Mw (this study) - Ms (GA)')
# plot 3: Log10M0 plot
# x = 3/2*Mw_ref+16.1
##x = [round(xx,2) for xx in Mw_ref]
##x = 1.5 * np.array(x) + 9.1
number1 = np.power(10,1.5*4.0 + 16.1)
number2 = np.power(10,1.5*3.9 + 16.1)
M0_ref = np.array([5.89e+21,1.50e+22,5.96e+22,5.37e+22,number1,4.57e+23,
                   1.12e+24,1.40e+22,3.13e+23,1.76e+24,4.07e+22,number2])
M0_ref = M0_ref * 1e-07
x = np.log10(M0_ref)
y = np.log10(np.array(M0_sim))
y_err = np.log10(np.array(M0_err_sim))
y_err = [y - y_err[:,0],y_err[:,1] - y]

fig = plt.figure(figsize=(10,10))
ax = plt.axes()
ax.plot(x,y,'o')
ax.errorbar(x,y,yerr = y_err,fmt='k.')
ax.plot([np.min(x)-0.5,np.max(x)+0.5],[np.min(x)-0.5,np.max(x)+0.5],'k-')
ax.set_xlabel('logMo (ref)')
ax.set_ylabel('logMo (this study)')

### plot 4: GOF leveles
##fig = plt.figure(figsize=(10,10))
##ax = plt.axes()
##xlabels = id_eqs
##ax.plot(range(len(misfit_sim)),misfit_sim,'ko')
##
##patches = []
##rect1 = mpatch.Rectangle((-0.5,0),len(xlabels),0.4)
##rect2 = mpatch.Rectangle((-0.5,0.4),len(xlabels),0.1)
##rect3 = mpatch.Rectangle((-0.5,0.5),len(xlabels),0.1)
##rect4 = mpatch.Rectangle((-0.5,0.6),len(xlabels),0.4)
##patches.append(rect1)
##patches.append(rect2)
##patches.append(rect3)
##patches.append(rect4)
##color = ['green','blue','yellow','red']
##collection = PatchCollection(patches,facecolors=color, alpha=0.3)
##ax.add_collection(collection)
##
##ax.set_xticks(range(len(misfit_sim)))
##ax.set_xticklabels(xlabels,fontsize=12)
##ax.set_xlim([-0.5,len(misfit_sim)-0.5])
plt.show()
    
