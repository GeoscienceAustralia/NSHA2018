# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to map

Usage:
    python map_nsha18.py <path to csv file>
    

@author: tallen
"""
from sys import argv
from matplotlib.mlab import griddata
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor
from tools.mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap # drawshapepoly, labelpolygon, 
import shapefile
from scipy.constants import g
from misc_tools import remove_last_cmap_colour
from shapely.geometry import Point, Polygon
import warnings
warnings.filterwarnings("ignore")
#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

drawshape = False # decides whether to overlay seismic sources

# set map resolution
res = 'i' 

# get current working directory (and computer!)
cwd = getcwd()

##############################################################################
# define inputs
##############################################################################

# set map file to plot
gridfile = argv[1]

# assign map name for plotting
modelName = argv[2]

# add contours?
addContours = argv[3] # True or False

# which probability - acceptable values are: 2 (2%), 9 (9.5%) or 10 (10%)
pltProbability = argv[4]

# plt GSHAP colours 
pltGSHAP = argv[5]  # True or False

##############################################################################
# parse hazard map file
##############################################################################

# parse hazard grid file 
lines = open(gridfile).readlines()

# get keys for model
if lines[0].startswith('#'):
    line = lines[1]
else:
    line = lines[0]

# get dictionary keys
keys = line.strip().split(',')[2:]

# make grid dictionary
grddict = []

print '\nReading', modelName
for line in lines[2:]:
    dat = line.strip().split(',')
    '''
    # check if GSHAP    
    if len(dat) == 1:
        dat = line.strip().split('\t')
        gshap = True
        keys = ['PGA-0.1']
    '''
    
    tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
    
    # fill keys
    idx = 2
    for key in keys:
        if pltGSHAP == 'True':
            # convert to m/s**2
            tmpdict[key] = float(dat[idx]) * g
        else:
            tmpdict[key] = float(dat[idx])
        idx += 1
    
    # add to grid list
    grddict.append(tmpdict)
    
##############################################################################    
# get key index for plotting
##############################################################################

for i, key in enumerate(keys):
    keyProb = str(int(floor(100*float(key.split('-')[-1]))))
    if keyProb == pltProbability:
        mapidx = i

##############################################################################    
# now make maps
##############################################################################

#keys = ['PGA_10', 'PGA_02', 'SA02_10', 'SA02_02', 'SA10_10', 'SA10_02']
for i, key in enumerate([keys[mapidx]]): # just plot 1 for now!
    if i > 0:
        plt.clf()
        plt.cla()

    # get IM period
    period = key.split('-')[0]
    period = period.replace('(','')
    period = period.replace(')','')
    period = period.replace('.','')
    
    print period
    
    # get map probability of exceedance
    probFraction = str(float(key.split('-')[-1]))
    probability = str(100*float(key.split('-')[-1])).split('.')[0]+'%'
    #probability = str(100*float(key.split('-')[-1]))+'%'
    if probability == '9%':
        probability = '9.5%'
    print 'Probability', probability
    
    figure = plt.figure(i,figsize=(19,12))
    
    ax = figure.add_subplot(111)
    
    bbox = '108/152/-44/-8' # map boundary - lon1/lon2/lat1/lat2
    bbox = '107.0/153.0/-45.0/-7.0'
    
    bbox = bbox.split('/')
    minlon = float(bbox[0])
    maxlon = float(bbox[1])
    minlat = float(bbox[2])
    maxlat = float(bbox[3])
    mbuff_l = 1.
    mbuff_r = 3.5
    
    '''
    # get shpfile for masking hazard values
    shpfile = solfile.split('.')[0]
    
    inshape = '../Grids/2005_grids/canada_2005grid_released.shp'
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    maskpoly = Polygon(sf[0].points)
    '''
    
    # build data to plot
    hazvals = []
    latlist = []
    lonlist = []
    
    # add buffer to data
    for gridval in grddict:
        lonlist.append(gridval['lon'])
        latlist.append(gridval['lat'])
        if gridval[key] == 0.0:
            hazvals.append(0.0)
        else:
            hazvals.append(gridval[key])
            
        '''
        # mask grid points outside defined grid to avoid extrapolation - this is slow!
        point = Point(gridval['lon'], gridval['lat'])
        if point.within(maskpoly) == False:
            hazvals.append(nan)
        else:
            hazvals.append(gridval[key])
        '''
    
    #idx = array(range(0, len(lonlist), 10)) # resample for quickly testing mapping
    idx = array(range(0, len(lonlist), 1))
    lonlist = array(lonlist)[idx]
    latlist = array(latlist)[idx]
    hazvals = array(hazvals)[idx]
    
    # delete zero hazvals
    idx =where(hazvals==0)[0]
    '''
    lonlist = delete(lonlist, idx)
    latlist = delete(latlist, idx)
    hazvals = delete(hazvals, idx)
    '''
    hazvals[idx] = 1E-20
    
    # get map bounds
    llcrnrlat = minlat
    urcrnrlat = maxlat
    llcrnrlon = minlon
    urcrnrlon = maxlon
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lon_0 = 134.
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    # set map
    # Projection used for National Mapping
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution=res,area_thresh=2000.)
                
    #m.drawmapboundary(fill_color='lightgray')
    #m.fillcontinents(color='white',lake_color='lightgray',zorder=0)
    m.drawcoastlines(linewidth=0.5,color='k')
    m.drawcountries(color='0.2')
    m.drawstates(color='0.2')
    
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
    extent = (minlon-mbuff_l, maxlon+mbuff_r, minlat-mbuff_r, maxlat+0)
    xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    	
    #resampled = griddata(lonlist, latlist, log10(hazvals), xs, ys, interp='linear')
    resampled = griddata(lonlist, latlist, hazvals, xs, ys, interp='linear')
    #resampled = griddata(lonlist, latlist, log10(hazvals), lonlist, latlist, interp='linear') # if this suddenly works, I have no idea why!
    
    # get 1D lats and lons for map transform
    lons = ogrid[extent[0]:extent[1]:N]
    lats = ogrid[extent[2]:extent[3]:N]
    
    # transform to map projection
    nx = int((m.xmax-m.xmin)/3000.)+1
    ny = int((m.ymax-m.ymin)/3000.)+1
    
    # differences in the way different machines deal with grids - weird!
    if cwd.startswith('/nas'):
        transhaz = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    else:
        transhaz = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    
    masked_array = ma.array(transhaz, mask=isnan(transhaz))
    
    # get colormap from cpt file
    cptfile = 'cw1-013_mod.cpt'
    #ncols = 9
    
    # get T from period
    if period.startswith('SA0'):
        T = 'Sa(0.'+period[3:]+')'
    elif period.startswith('SA'):
        T = 'Sa('+period[2]+'.'+period[3:]+')' 
    else:
        T = period
        
    print period, T
    
    if probability == '10%' or probability == '9.5%': # kluge to get on same scale
        ncolours = 13
        
    elif probability == '2%':
        ncolours = 13       
    
    ncolours = 13
    try:
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=True)
    except:
        try:
            if pltGSHAP == 'True':
                nascptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gshap_mpl.cpt'
                capfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/capitals_names.csv'
                ncolours = 10
                cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=False)
                cmap = remove_last_cmap_colour(cmap)
                
            else:
                nascptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/'+ cptfile
                capfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/capitals_names.csv'
                cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=True)
                #cmap = remove_last_cmap_colour(cmap)
            
            #cptfile = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/GMT_no_green.cpt'
            
        except:
            try:
                ncicptfile = '/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/maps/'+ cptfile
                capfile = '/short/w84/NSHA18/sandpit/tia547/NSHA2018/shared/capitals_names.csv'
                cmap, zvals = cpt2colormap(ncicptfile, ncolours, rev=True)

            except:
                ncicptfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/maps/'+ cptfile
                capfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/capitals_names.csv'
                cmap, zvals = cpt2colormap(ncicptfile, ncolours, rev=True)
    
    print 'Making map...'    
    cmap.set_bad('w', 1.0)
    
    if probability == '10%' or probability == '9.5%':
        if pltGSHAP == 'True':
            bounds = array([0., 0.2, 0.4, 0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 6.0])
            #ncolours = 9
            #norm = colors.Normalize(vmin=0,vmax=10)
        else:
            if period == 'PGA' or period == 'SA005' or period == 'SA01' or period == 'SA02' \
               or period == 'SA03' or period == 'SA05':
                bounds = array([0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
            elif  period == 'SA07'  or period == 'SA10':
                bounds = array([0, 0.002, 0.004, 0.007, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08])
            else:
                bounds = array([0, 0.001, 0.002, 0.004, 0.007, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.045, 0.06])
            ncolours = 12
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    else:
        if period == 'PGA' or period == 'SA005' or period == 'SA01' or period == 'SA02' \
           or period == 'SA03' or period == 'SA05':
            bounds = array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.3, 0.5, 0.7, 1.0])
        elif period == 'SA07':
            bounds = array([0, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24, 0.36])
        elif period == 'SA15' or period == 'SA10':
            bounds = array([0, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.2, 0.3])
        else:
            bounds = array([0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.015, 0.03, 0.04, 0.05, 0.06, 0.1, 0.16])
        #ncolours = 12
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    
    #m.imshow(masked_array, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax, zorder=0)
    m.imshow(masked_array, cmap=cmap, extent=extent, norm=norm, zorder=0)
    
    ##########################################################################################
    # plot contours
    ##########################################################################################
    if addContours == 'True':
        x, y = m(xs, ys)
        if probability == '10%':
            '''
            levels = arange(0.02, 0.3, 0.02)
            levels = arange(0.05, 0.3, 0.05)
            levels = array([0.01, 0.02, 0.04, 0.06, 0.08001, 0.12, 0.16, 0.18, 0.24])
            '''
            levels_lo = bounds[1:4]
            levels = bounds[4:]
        elif probability == '2%':
            #levels = arange(0.05, 0.31, 0.05)
            levels = bounds[1:]
        
        csm = plt.contour(x, y, resampled, levels, colors='0.2', lw=0.3)    
        csm_lo = plt.contour(x, y, resampled, levels_lo, colors='0.2', lw=0.3)
        
        plt.clabel(csm, inline=1, fontsize=10, fmt='%0.2f')
        plt.clabel(csm_lo, inline=1, fontsize=10, fmt='%0.3f')
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    # mask non-AU polygons
    nonmask = [0, 2, 3, 4, 6, 7, 11, 13, 16, 17] # polygon number
    landpolys = []
    for pidx, polygon in enumerate(m.landpolygons):
        maskPoly = True
        for nmidx in nonmask:
            if pidx == nmidx:
                maskPoly = False 
        if maskPoly == True:
            poly = polygon.get_coords()
            plt.fill(poly[:,0], poly[:,1], 'w')
        
    #mask_outside_polygon(polys[1][::-1], ax=None)
    polys = get_map_polygons(m)
    mask_outside_polygons(polys, '0.9', plt)
    
    # get lake ploygons
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], '0.9')
        polygons.append(poly)
    
    ##########################################################################################
    # format main axis
    ##########################################################################################
    if probability == '9.5%':
        titlestr = ' '.join((modelName, '1 in 500-Year AEP Mean', T, 'Hazard on AS1170.4 Site Class '))
    else:
        titlestr = ' '.join((modelName, T, probability, 'in 50-Year Mean Hazard on AS1170.4 Site Class '))    
    
    # comment out for final GA Record
    #plt.title(titlestr+'$\mathregular{B_e}$')
    
    # get map bbox
    if i == 0:
        map_bbox = ax.get_position().extents
    
    ##########################################################################################
    # add DRAFT text!
    ##########################################################################################
    '''
    import matplotlib.patheffects as path_effects
    #import matplotlib.patheffects as PathEffects
    drafttext = figure.text(0.5, 0.5, 'DRAFT', color='w', rotation=45,
                          ha='center', va='center', size=160, alpha=0.1)
    #drafttext = ax.annotate("DRAFT", xy=(.5, .5), xytext=(.5, .5),
    #                        ha='center', va='center', size=160, rotation=45)                         
    drafttext.set_path_effects([path_effects.Stroke(linewidth=4, foreground='maroon'),
                       path_effects.Normal()])
    drafttext.set_alpha(0.5)
    ''' 
    
    '''
    ###########################################################################################
    annotate cities
    ###########################################################################################
    '''
    capfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/capitals_names.csv'
    
    import matplotlib.patheffects as PathEffects
    pe = [PathEffects.withStroke(linewidth=2.5, foreground="w")]
              
    llat = []
    llon = []
    locs = []
    textoffset = []
    
    # read data
    lines = open(capfile).readlines()
    for line in lines:
        llon.append(float(line.strip().split(',')[0]))
        llat.append(float(line.strip().split(',')[1]))
        locs.append(line.strip().split(',')[2])
        textoffset.append(float(line.strip().split(',')[3]))
    
    # plot locs on map
    x, y = m(array(llon), array(llat))
    plt.plot(x, y, 's', markerfacecolor='None', markeredgecolor='k', markeredgewidth=0.5, markersize=8)
    
    # label cities
    for i, loc in enumerate(locs):
        if textoffset[i] == 0.:
            x, y = m(llon[i]-0.35, llat[i]+0.12)
            plt.text(x, y, loc, size=15, ha='right', va='bottom', weight='light')
        elif textoffset[i] == 1.:
            x, y = m(llon[i]+0.35, llat[i]+0.12)
            #plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
            plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light')
        elif textoffset[i] == 2.:
            x, y = m(llon[i]+0.3, llat[i]-0.3)
            plt.text(x, y, loc, size=15, ha='left', va='top', weight='light')
        elif textoffset[i] == 3.:
            x, y = m(llon[i]-0.3, llat[i]-0.2)
            plt.text(x, y, loc, size=15, ha='right', va='top', weight='light')
    
    ##########################################################################################
    # add GA logo
    ##########################################################################################
    if modelName.startswith('GSHAP') == False:
        # load logo
        try:
            im = plt.imread('../GAlogo.png')
        except:
            # cover all bases
            try:
                im = plt.imread('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/GAlogo.png')
            except:
                try:
                    im = plt.imread('/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/GAlogo.png')
                except:
                    im = plt.imread('/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/GAlogo.png')
        
        # set bbox for logo
        imoff = 0.02
        #logo_bbox = mpl.transforms.Bbox(array([[map_bbox[0]+imoff,map_bbox[1]+imoff],[0.15,0.15]]))
        #logo_bbox = [map_bbox[0]+0.11,map_bbox[1]-0.005,0.15,0.15]
        logo_bbox = [map_bbox[0]+0.14,map_bbox[1]-0.075,0.25,0.25]
        newax = figure.add_axes(logo_bbox) #, zorder=-1)
        newax.imshow(im)
        newax.axis('off')
        
        ##########################################################################################
        # add CC-by
        ##########################################################################################
        
        # load logo
        try:
            im = plt.imread('../ccby_narrow.png')
        except:
            # covering all bases again
            try:
                im = plt.imread('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/ccby_narrow.png')
                
            except:
                try:
                    im = plt.imread('/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/ccby_narrow.png')
                    
                except:
                    im = plt.imread('/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/ccby_narrow.png')
                    
    
        # set bbox for logo
        imoff = 0.02
        logo_bbox = [map_bbox[0]+0.11,map_bbox[1]-0.005,0.2,0.2]
        logo_bbox = [0.66,map_bbox[1]-0.03,0.1,0.1]
        newax = figure.add_axes(logo_bbox) #, zorder=-1)
        newax.imshow(im)
        newax.axis('off')
    
        
    
    ##########################################################################################
    # superimpose area source shapefile
    ##########################################################################################
    '''
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
    make colourbar
    ###########################################################################################
    '''    
    
    # set colourbar
    plt.gcf().subplots_adjust(bottom=0.1)
    cax = figure.add_axes([0.31,0.05,0.38,0.02]) # setup colorbar axes.
    #norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
    
    # set cb labels
    #linticks = array([0.01, 0.03, 0.1, 0.3 ])
    #logticks = arange(vmin, vmax+0.25, 0.25)
    #cb.set_ticks(logticks)
    #labels = [str('%0.3f' % 10**x) for x in logticks]
    
    cb.set_ticks(bounds)
    if pltGSHAP == 'True':
        labels = [str('%0.1f' % x) for x in bounds]
    else:
        labels = ['0'+str('%0.3f' % x).strip('0') for x in bounds]
    labels[0] = '0.0'
    if labels[-1] == '01.':
        labels[-1] = '1.0'
    
    cb.set_ticklabels(labels)
    
    # set title
    if pltGSHAP == 'True':
        titlestr = 'PGA 10% in 50-Year Mean Hazard ($\mathregular{m/s_2}$)'
    else:
        if probability == '9.5%':
            titlestr = ' '.join(('1 in 500-Year AEP Mean', T, 'Hazard (g)'))
        else:
            titlestr = ' '.join((T, probability, 'in 50-Year Mean Hazard (g)'))
    cb.set_label(titlestr, fontsize=12)
    
    # check to see if maps exists
    if path.isdir('maps') == False:
        mkdir('maps')
        
    # now save png file
    plt.savefig(path.join('maps', 'hazard_map_'+modelName.replace(' ','_')+'.'+period+'.'+probFraction+'.png'), \
                dpi=300, format='png', bbox_inches='tight')
    
    # save pdf file
    '''
    plt.savefig(path.join('maps', 'hazard_map_'+modelName.replace(' ','_')+'.'+key+'.pdf'), \
                dpi=300, format='pdf', bbox_inches='tight')
    '''
    #plt.show()
    plt.close('all')
    
    ##########################################################################################
    # make shapefile of contour lines
    ##########################################################################################
    
    print 'Masking maritime boundaries...'
    if getcwd().startswith('/nas'):
        inshape = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/shapefiles//au_maritime_boundary_digitised.shp'
    else:
        inshape = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/maps/shapefiles/au_maritime_boundary_digitised.shp'
    
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    poly = Polygon(sf[0].points)
    
    
    # check to see if shapefile contours exists
    if path.isdir('contours') == False:
        mkdir('contours')
        
    # make list of levels - old levels array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.12, 0.18, 0.24])
    allLevels = [bounds] #,
    '''
                 array([0.03, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.25, 0.30]),
                 array([0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]),
                 array([0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.15]),
                 array([0.003, 0.004, 0.006, 0.008, 0.01, 0.012, 0.015, 0.02])]
    '''
    levelNames = ['lev_nat', 'lev_swsz','lev_ntsa', 'lev_nswtas', 'lev_qld'] #, 'lev_0_01', 'lev_0_02', 'lev_0_05']  
    levelNames = ['lev_nat']
    
    # loop thru levels
    for levels, levelName in zip(allLevels, levelNames):
        
        # setup shapefile
        outshp = path.join('contours', '_'.join((modelName.replace(' ','_'), key, \
                           levelName, 'contours.shp')))
    
        # set shapefile to write to
        w = shapefile.Writer(shapefile.POLYLINE)
        w.field('LEVELS','F', 5, 3)
            
        # have to re-contour using un-transformed lat/lons
        cs = plt.contour(xs, ys, resampled, levels, colors='k')
        
        plt.close(figure)
        
        # loop through contour levels
        for l, lev in enumerate(cs.levels):
            contours = cs.collections[l].get_paths()
            
            # now loop through multiple paths within level
            for cnt in contours:
                # check if verticies within maritime boundaries
                newcnt = []
                for vert in cnt.vertices:
                    point = Point(vert)
                    if point.within(poly) == True:
                        newcnt.append(vert)
                    
                    # add polyline to shapefile
                    elif point.within(poly) == False and len(newcnt) > 0:
                        #w.line(parts=[newcnt.vertices], shapeType=shapefile.POLYLINE)
                        w.line(parts=array([newcnt]), shapeType=shapefile.POLYLINE)                        
                        # add level attribute
                        w.record(lev)
                        
                        newcnt = []
                        
                # do final addition of contours
                if len(newcnt) > 0:
                    w.line(parts=array([newcnt]), shapeType=shapefile.POLYLINE)
                    # add level attribute
                    w.record(lev)
                
        # now save area shapefile
        w.save(outshp)
        
        # write projection file
        prjfile = outshp.strip().split('.shp')[0]+'.prj'
        f = open(prjfile, 'wb')
        f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
        f.close()
        

    


















