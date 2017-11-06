from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

cat_file = './catalogue.txt'
data_dir = './output/'

plt.figure('1')
# draw the background
m = Basemap(projection='cyl',resolution=None,llcrnrlon=110.,
            llcrnrlat=-45.,urcrnrlon=160.,urcrnrlat= -10.)

m.etopo()

cat = np.loadtxt(cat_file)


for e in cat:
    station_file = data_dir + 'stations_eve_' + str(int(e[0])) + '.txt'
    with open(station_file) as f:
        lines = f.readlines()
    x_eve,y_eve = m(e[2],e[1])
    for line in lines:
        line = line.split(' ')
        x_sta,y_sta = m(float(line[2]),float(line[1]))
        dist_sta = float(line[3]) * 111.194929703
        plt.figure('1')
        m.scatter(x_sta,y_sta,50,marker='^',color='k')
        m.plot([x_eve,x_sta],[y_eve,y_sta],color='#999999',alpha=0.99)
        
        plt.figure('2')
        plt.scatter(dist_sta,e[13],80,marker='o',color='b',alpha=0.5)
        plt.xlim([0,1200])
        plt.ylim([3,6])
        plt.xlabel('Distance(km)')
        plt.ylabel('Magnitude(Ml)')
        
plt.figure('1')        
x_eve,y_eve = m(cat[:,2],cat[:,1])    
m.scatter(x_eve,y_eve,80,marker='h',color='r')

##x_b,y_b = m(116.7091,-30.6147 )
##m.scatter(x_b,y_b,100,marker='o',color='b')






plt.show()
