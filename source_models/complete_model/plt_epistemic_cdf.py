import matplotlib.pyplot as plt
from numpy import arange, array, interp, argsort
from os import path, getcwd, sep
from sys import argv
import matplotlib as mpl
from misc_tools import dictlist2array
mpl.style.use('classic')

fracFolder = argv[1] # folder where fractile files sit

altPlaces = True
    
# set params 
fractiles = arange(0., 1.01, 0.01)

# get keys to fill
fracFile = path.join(fracFolder, 'quantile_map-0.5_1.csv' )

# parse file
lines = open(fracFile).readlines()

# get keys from medianfile
keys = lines[1].strip().split(',')[2:]
	
pltlett = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']

	
# load mean lines to add
meanFile = path.join(fracFolder, 'hazard_map-mean_1.csv' )
meanLines = open(meanFile).readlines()[2:]

# loop thru fractiles and fill city Dict
fracDict = []
for i, fractile in enumerate(fractiles):
    fracFile = path.join(fracFolder, 'quantile_map-'+str(round(fractile,2))+'_1.csv' )
    
    # parse file
    lines = open(fracFile).readlines()[2:]
    
    # loop through lines
    j = 0
    for line, mLine in zip(lines, meanLines):
        dat = line.strip().split(',')
        mdat = mLine.strip().split(',')
        
        # if first fractile
        if i == 0:
            tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
            for k, key in enumerate(keys):
                tmpdict['quant_'+key] = [float(dat[k+2])]
                
                # add mean values too
                tmpdict['mean_'+key]  = float(mdat[k+2])

            fracDict.append(tmpdict)
            
        # now append to existing dict
        else:
            for k, key in enumerate(keys):
                fracDict[j]['quant_'+key].append(float(dat[k+2]))
                
        j+=1
                
# now turn lists to arrays
for j in range(0, len(fracDict)):
    for k, key in enumerate(keys):
        fracDict[j]['quant_'+key] = array(fracDict[j]['quant_'+key])
        
###################################################################################
# match city name to fracDict
###################################################################################

# first parse city file
cwd = getcwd()
if cwd.startswith('/Users'): #mac
    citycsv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
else:
    citycsv = '../../shared/nsha_cities.csv'
lines = open(citycsv).readlines()
    
# make city dict
cityDict = []
for line in lines:
    dat = line.strip().split(',')
    tmpdict = {'city':dat[2], 'lon':float(dat[0]), 'lat':float(dat[1])} 
    cityDict.append(tmpdict)


for j in range(0, len(fracDict)):
    for city in cityDict:
        if city['lon'] == fracDict[j]['lon'] \
           and city['lat'] == fracDict[j]['lat']:
           
           # add place
           fracDict[j]['place'] = city['city']
           
###################################################################################
# export csvs
###################################################################################
fracPath = path.join('final', 'fractile_tables')

# get alphabetical order
placesOrd = dictlist2array(fracDict,'place')
sortidx = argsort(placesOrd)
header = 'PLACE,LONGITUDE,LATITUDE,MEAN,16TH PERCENTILE,50TH PERCENTILE,84TH PERCENTILE,MEAN PERCENTILE\n'

for key in keys:
    csvfile = 'fractile_table_' + key + '.csv'
    fracFilePath = path.join(fracPath, csvfile)
    
    pctltxt = header
    
    for idx in sortidx:
        # get fractile at mean
        if not fracDict[idx]['place'].endswith('max'):
            meanfrac = interp(fracDict[idx]['mean_'+key], fracDict[idx]['quant_'+key], fractiles)
            
            pctltxt = pctltxt+','.join((fracDict[idx]['place'], \
                                str('%0.2f' % fracDict[idx]['lon']), \
                                str('%0.2f' % fracDict[idx]['lat']), \
                                str('%0.3f' % fracDict[idx]['mean_'+key]), \
                                str('%0.3f' % fracDict[idx]['quant_'+key][16]), \
                                str('%0.3f' % fracDict[idx]['quant_'+key][50]), \
                                str('%0.3f' % fracDict[idx]['quant_'+key][84]), \
                                str(int(round(100*meanfrac))))) + '\n'
        
    # write to file
    f = open(fracFilePath, 'wb')
    f.write(pctltxt)
    f.close()

###################################################################################
# let's make the plots
###################################################################################
if altPlaces == False:
    places = ['Perth', 'Darwin', 'Adelaide', 'Melbourne', 'Hobart', 'Canberra', 'Sydney', 'Brisbane']
else:
    places = ['Wongan Hills', 'Kalgoorlie', 'Port Pirie', 'Yulara', 'Tennant Creek', 'Cooma', 'Albury', 'Morwell']

# loop through keys
xtxt = 0.28
ytxt = 0.97
# just do PGA
for k, key in enumerate(keys[:3]):
    # set up figure
    fig = plt.figure(k+1, figsize=(15, 15))
    
    # loop thru places to plot
    for i, place in enumerate(places):
        # loop thru fracDict
        for frac in fracDict:
            if place == frac['place']:
            
                # plot fig
                ax = plt.subplot(4, 2, i+1)
                plt.semilogx(frac['quant_'+key], fractiles, 'k-', lw=1.5, label='Fractiles')
                
                
                # make pretty
                plt.title(place, fontsize=18)
                if i == 0 or i == 2 or i == 4 or i == 6:
                    plt.ylabel('Fractile', fontsize=16)
                
                if i >= 6:
                    plt.xlabel(key.replace('(','').replace(')','').split('-')[0] + ' (g)', fontsize=16)
                
                plt.grid(which='both')
                
                if key.endswith('0.02'):
                    plt.xlim([0.01, 0.8])
                else:
                    plt.xlim([0.003, 0.3])
                
                # plt mean
                plt.semilogx([frac['mean_'+key],frac['mean_'+key]], [0,1], '--', c='seagreen', lw=1.5, label='Mean')
                # plt median
                plt.semilogx([frac['quant_'+key][50],frac['quant_'+key][50]], [0,1], '--', c='dodgerblue', lw=1.5, label='50th Percentile')
                # plt 84th
                plt.semilogx([frac['quant_'+key][84],frac['quant_'+key][84]], [0,1], '--', c='orange', lw=1.5, label='84th Percentile')
                # plt 95th
                plt.semilogx([frac['quant_'+key][95],frac['quant_'+key][95]], [0,1], '--', c='r', lw=1.5, label='95th Percentile')
                
                if i == 0:
                    plt.legend(loc=2, fontsize=13)
                    
                plt.text(xtxt, ytxt, pltlett[i], fontsize=18, va='top', ha='right')
                
    #plt.suptitle(fracFolder.split(sep)[1] + ' ' + key, fontsize=20)
    
    # set fig file
    if altPlaces == True:
        figFile = path.join(fracFolder, '_'.join((fracFolder.split(sep)[1],key,'alt_CDF.png')))
    else:
        figFile = path.join(fracFolder, '_'.join((fracFolder.split(sep)[1],key,'CDF.png')))
    plt.savefig(figFile, fmt='png', bbox_inches='tight')
plt.show()