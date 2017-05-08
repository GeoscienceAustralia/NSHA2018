# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to shp

Usage:
    python sol2map.py minlon/maxlon/minlat/maxlat
    
    eg:
    
    python sol2map.py -140/-114/52/70

@author: tallen
"""
from sys import argv
from matplotlib.mlab import griddata
from matplotlib import colors, colorbar, cm
from os import path, mkdir
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan
from tools.mapping_tools import drawshapepoly, labelpolygon, get_map_polygons, mask_outside_polygons, cpt2colormap

import shapefile
#from gmt_tools import cpt2colormap
from shapely.geometry import Point, Polygon

mpl.rcParams['pdf.fonttype'] = 42

drawshape = True
#bbox = '-133/-120/48/56'

# set map file to plot
gridfile = argv[1]

# get model
model = path.split(gridfile)[-1].split('_')[2] # this will likely need modifying depending on filename format

# set map resolution
res = 'l' 

# parse sol file 
lines = open(gridfile).readlines()

# make grid dictionary
grddict = []

# get keys for model
if lines[0].startswith('#'):
    line = lines[1]
else:
    line = lines[0]

# get dictionary keys
keys = line.strip().split(',')[2:]

print '\nReading', model
for line in lines[2:]:
    tmpdict = {}
    dat = line.strip().split(',')
    tmpdict['lon'] = float(dat[0])
    tmpdict['lat'] = float(dat[1])
    
    # fill keys
    idx = 2
    for key in keys:
        tmpdict[key] = float(dat[idx])
        idx += 1
    
    # add to grid list
    grddict.append(tmpdict)
    
'''    
# now make maps
'''
#keys = ['PGA_10', 'PGA_02', 'SA02_10', 'SA02_02', 'SA10_10', 'SA10_02']

for i, key in enumerate([keys[1]]): # just plot 1!
    
    # get IM period
    period = key.split('-')[0]
    
    # get map probability of exceedance
    probability = str(100*float(key.split('-')[-1])).split('.')[0]+'%'
    
    #figure = plt.figure(i,figsize=(19,12))
    plt.clf()
    plt.cla()
    figure, ax = plt.subplots(i+1,figsize=(19,12))
    #ax = plt.subplots(111)
    #bbox = argv[1].split('/')
    
    bbox = '108/153/-44/-8'
    bbox = bbox.split('/')
    minlon = float(bbox[0])
    maxlon = float(bbox[1])
    minlat = float(bbox[2])
    maxlat = float(bbox[3])
    mbuff = 1.
    
    # build data to plot
    hazvals = []
    latlist = []
    lonlist = []
    
    '''
    # get shpfile name
    shpfile = solfile.split('.')[0]
    
    inshape = '../Grids/2005_grids/canada_2005grid_released.shp'
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    poly = Polygon(sf[0].points)
    '''
    
    # add buffer to data
    for gridval in grddict:
        lonlist.append(gridval['lon'])
        latlist.append(gridval['lat'])
        if gridval[key] == 0.0:
            hazvals.append(0.0)
        else:
            hazvals.append(gridval[key])
            
        '''
        # mask grid points outside defined grid to avoid extrapolation
        point = Point(gridval['lon'], gridval['lat'])
        if point.within(poly) == False:
            hazvals.append(nan)
        else:
            hazvals.append(gridval[key])
        '''
    
    idx = array(range(0, len(lonlist), 100))
    idx = array(range(0, len(lonlist), 1))
    lonlist = array(lonlist)[idx]
    latlist = array(latlist)[idx]
    hazvals = array(hazvals)[idx]
    
    # get map bounds
    llcrnrlat = minlat
    urcrnrlat = maxlat
    llcrnrlon = minlon
    urcrnrlon = maxlon
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    # set map
    # Projection used for National Mapping
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution='l',area_thresh=1000.)
                
    #m.drawmapboundary(fill_color='lightgray')
    #m.fillcontinents(color='white',lake_color='lightgray',zorder=0)
    m.drawcoastlines(linewidth=0.5,color='k')
    m.drawcountries(color='0.2')
    m.drawstates(color='0.5')
    
    # draw parallels and meridians.
    if maxlon-minlon > 40:
        xlabel = 6.
    elif maxlon-minlon > 20:
        xlabel = 4.
    elif maxlon-minlon > 10:
        xlabel = 2.
    else:
        xlabel = 1.
        
    if maxlat-minlat > 40:
        ylabel = 6.
    elif maxlat-minlat > 20:
        ylabel = 4.
    elif maxlat-minlat > 10:
        ylabel = 2.
    else:
        ylabel = 1.
            
    m.drawparallels(arange(-90.,90.,ylabel), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    m.drawmeridians(arange(0.,360.,xlabel), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    
    # first make regular cartesian grid
    print 'Resampling data...'
    N = 500j
    extent = (minlon-mbuff, maxlon+mbuff, minlat-mbuff, maxlat+mbuff)
    xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    resampled = griddata(lonlist, latlist, log10(hazvals), xs, ys, interp='linear')
    
    
    # get 1D lats and lons for map transform
    lons = ogrid[extent[0]:extent[1]:N]
    lats = ogrid[extent[2]:extent[3]:N]
    
    # transform to map projection
    nx = int((m.xmax-m.xmin)/2000.)+1
    ny = int((m.ymax-m.ymin)/2000.)+1
    transhaz = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    masked_array = ma.array(transhaz, mask=isnan(transhaz))
    #masked_array = masked_array.set_fill_value(0)
    
    # get colormap from cpt file
    cptfile = 'cw1-013.cpt'
    ncols = 9
    
    #cmap = cm.rainbow
    if period == 'PGA':
        
        if probability == '10%':
            ncolours = 14
            vmin = -2.
            vmax = -0.25
            
        elif probability == '2%':
            ncolours = 12
            vmin = -1.5
            vmax = vmin + 0.25 * ncolours/2.
        T = 'PGA'
        
    elif period == 'SA02':
        ncolours = 14
        if probability == '10%':
            vmin = -3
            vmax = vmin + 0.5 * ncolours/2.
        
        elif probability == '2%':
            vmin = -1.75
            vmax = vmin + 0.25 * ncolours/2.
        T = 'Sa(0.2 s)'
        
    elif period == 'SA10':
        
        if probability == '10%':
            vmin = -3
            vmax = vmin + 0.25 * ncolours/2.
        
        elif probability == '2%':
            ncolours = 14
            vmin = -2
            vmax = vmin + 0.25 * ncolours/2.
        T = 'Sa(1.0 s)'
    
    try:
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=True)
    except:
        cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/precip3_16lev.cpt'
        #cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
    
    print 'Making map...'    
    cmap.set_bad('w', 1.0)
    m.imshow(masked_array, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax, zorder=0)
    
    ##########################################################################################
    # plot contour
    ##########################################################################################
    
    x, y = m(xs, ys)
    if probability == '10%':
        levels = arange(0.02, 0.3, 0.02)
    elif probability == '2%':
        levels = arange(0.05, 0.3, 0.05)
    CS = m.contour(x, y, 10**resampled, levels, colors='k')
    
    plt.clabel(CS, inline=1, fontsize=10)
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    polys = get_map_polygons(m)
    
    #mask_outside_polygon(polys[1][::-1], ax=None)
    mask_outside_polygons(polys, '0.9', plt)
    
    # get lake ploygons
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
        polygons.append(poly)
    
    plt.title(' '.join((model, T, probability, 'in 50-Year Mean Hazard on Site Class B/C')))
    
    # get map bbox
    map_bbox = ax.get_position().extents
    
    ##########################################################################################
    # add GA logo
    ##########################################################################################
    
    # load logo
    im = plt.imread('../logo.png')
    
    # set bbox for logo
    imoff = 0.02
    logo_bbox = mpl.transforms.Bbox(array([[map_bbox[0]+imoff,map_bbox[1]+imoff],[0.15,0.15]]))
    logo_bbox = [map_bbox[0]+0.11,map_bbox[1]-0.005,0.15,0.15]
    newax = figure.add_axes(logo_bbox) #, zorder=-1)
    newax.imshow(im)
    newax.axis('off')
     
    
    '''
    # superimpose area source shapefile
    shpfile =['..//final_inputs//SWCan_T3EclC_area1.shp', \
              '..//final_inputs//WArctic_area.shp']
    if drawshape == True:
        for shp in shpfile:
            sf = shapefile.Reader(shp)
            drawshapepoly(m, plt, sf, col='k', lw=1.5, polyline=True)
            labelpolygon(m, plt, sf, 'CODE', fsize=14)    
    '''
    '''
    ###########################################################################################
    annotate probabilities
    ###########################################################################################
    '''
    """
    import matplotlib.patheffects as PathEffects
    pe = [PathEffects.withStroke(linewidth=2.5, foreground="w")]
    
    # set cities
    locs = ['Yellowknife', 'Inuvik', 'Whitehorse', 'Sandspit', 'Victoria', \
            'Vancouver', 'Kamloops', 'Cranbrook', 'Calgary', 'Edmonton', 'Regina', \
            'Winnipeg', 'Thunder Bay', 'Windsor', 'Sudbury', 'Toronto', 'Ottawa', \
            'Montreal', 'Quebec City', 'Riviere-du-Loup', 'Fredericton', 'Charlottetown', \
            'Halifax', "St. John's", 'Iqaluit']
            
    locidx = [621, 616, 606, 62, 88, 82, 30, 15, 96, 106, 167, 199, 393, 416, 387, \
              646, 338, 678, 497, 502, 543, 578, 561, 595, 632]
              
    tbflag = [1, 1, 1, 1, -2, 0, -3, 0, 1, 1, 1, 1, -2, -3, 1, 0, -2, 0, -3, -3, 0, 1, 1, -3, 1]
              
    pfile = 'MMI8.0_50-yr.locs.csv'
    
    llat = []
    llon = []
    
    # read data
    lines = open(pfile).readlines()[1:]
    for idx in locidx:
        llon.append(float(lines[idx].strip().split(',')[1]))
        llat.append(float(lines[idx].strip().split(',')[2]))
    
    # plot locs on map
    x, y = m(llon, llat)
    plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2., markersize=10)
    
    # label cities
    for i, loc in enumerate(locs):
        if tbflag[i] == 1:
            x, y = m(llon[i]+0.3, llat[i]+0.3)
            plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='normal', path_effects=pe)
        elif tbflag[i] == 0 or tbflag[i] == -1:
            x, y = m(llon[i]+0.3, llat[i]-0.35)
            plt.text(x, y, loc, size=15, ha='left', va='top', weight='normal', path_effects=pe)
        elif tbflag[i] == -2:
            x, y = m(llon[i]-0.3, llat[i]-0.35)
            plt.text(x, y, loc, size=15, ha='right', va='top', weight='normal', path_effects=pe)
        elif tbflag[i] == -3:
            x, y = m(llon[i]-0.35, llat[i]+0.3)
            plt.text(x, y, loc, size=15, ha='right', va='bottom', weight='normal', path_effects=pe)
    """
    
    '''
    ###########################################################################################
    make colourbar
    ###########################################################################################
    '''
    
    
    # set colourbar
    plt.gcf().subplots_adjust(bottom=0.1)
    cax = figure.add_axes([0.34,0.05,0.33,0.02]) # setup colorbar axes.
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
    
    # set cb labels
    #linticks = array([0.01, 0.03, 0.1, 0.3 ])
    logticks = arange(vmin, vmax+0.25, 0.25)
    cb.set_ticks(logticks)
    labels = [str('%0.2f' % 10**x) for x in logticks]
    cb.set_ticklabels(labels)
    
    # get map probabiltiy from filename
    probstr = '2'# yrs
    itime = '2,475'
    
    titlestr = ' '.join((T, probability, 'in 50-Year Mean Hazard (g)'))
    cb.set_label(titlestr, fontsize=12)
    
    # check to see if exists
    if path.isdir(key) == False:
        mkdir(key)
    
    plt.savefig(gridfile.strip('csv')+key+'.png', dpi=150, format='png', bbox_inches='tight')
    
plt.show()
