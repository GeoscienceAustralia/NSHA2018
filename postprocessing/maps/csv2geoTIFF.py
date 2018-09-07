# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:07:46 2015

@author: tallen
"""
'''
ncols=174
nrows=115
xllcorner=14.97
yllcorner=-34.54
cellsize=0.11
'''

'''
# read sol file
'''
from sys import argv
from os import sep
from numpy import array, mgrid, nan, shape
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
from osgeo import osr, gdal
import shapefile	
from shapely.geometry import Point, Polygon

period = argv[1]  # fmt: e.g. Sa1.0 or PGA

#solfile = ''.join(('solfiles', sep, 'NBCC2015_Canada_', period, 'all_000thperc.sol'))
solfile = ''.join(('solfiles', sep, 'GSC2015_seismichazardgrid_', period, '.txt'))

# parse sol file 
lines = open(solfile).readlines()[2:]

# make grid dictionary
grddict = []
alon = []
alat = []


maxlat = -90
minlat = 90
maxlon = -180
minlon = 180

print '\nReading sol file...'
for line in lines:
    tmpdict = {}
    dat = line.strip().split()
    tmpdict['lat'] = float(dat[0])
    tmpdict['lon'] = float(dat[1])
    tmpdict['P0.0200'] = float(dat[2])
    tmpdict['P0.01375'] = float(dat[3])
    tmpdict['P0.0100'] = float(dat[4])
    tmpdict['P0.00445'] = float(dat[5])
    tmpdict['P0.0021'] = float(dat[6])
    tmpdict['P0.0010'] = float(dat[7])
    tmpdict['P0.0005'] = float(dat[8])
    tmpdict['P0.000404'] = float(dat[9])
    #tmpdict['P0.0002'] = float(dat[10])
    #tmpdict['P0.0001'] = float(dat[11])
    tmpdict['ref'] = int(dat[11])
    
    grddict.append(tmpdict)
    alon.append(tmpdict['lon'])
    alat.append(tmpdict['lat'])
    
    # get bbox 
    if tmpdict['lon'] > maxlon:
        maxlon = tmpdict['lon']
    if tmpdict['lon'] < minlon:
        minlon = tmpdict['lon']
        
    if tmpdict['lat'] > maxlat:
        maxlat = tmpdict['lat']
    if tmpdict['lat'] < minlat:
        minlat = tmpdict['lat']

'''
make mesh
'''
resolution = 0.05 # degrees
invres = int(1/resolution)
keys = ['P0.0021', 'P0.0010', 'P0.000404'] # probabilities
for key in keys:
    print 'Making', key, 'grid mesh...'
    
    # first make z data array
    ahaz = []
    for grdval in grddict:
        ahaz.append(grdval[key])
    
    # space at define resolution    
    nlon = int(round((maxlon-minlon)*invres))
    nlat = int(round((maxlat-minlat)*invres))
    grid_x, grid_y = mgrid[minlon:maxlon:complex(nlon), minlat:maxlat:complex(nlat)]  # lon: 1000 @ 0.1 deg; lat: 500 @ 0.1 deg
    
    # interpolate z data array
    grid_z = griddata((array(alon), array(alat)), array(ahaz), (grid_x, grid_y), \
                      method='cubic', fill_value=nan) # scipy.interpolate (linear, nearest, cubic)
    #grid_z = griddata(array(alon), array(alat), array(ahaz), grid_x, grid_y, interp='linear') # matplotlib
    
    # mask grid points outside defined grid to avoid extrapolation
    print 'Masking', key, 'grid...'
    #inshape = '2015NBCC_grid_mask.shp'
    inshape = '..\\2005_grid\\canada_2005grid_released.shp'
    
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    poly = Polygon(sf[0].points)
    flat_x = grid_x.flatten()
    flat_y = grid_y.flatten()
    flat_z = grid_z.flatten()
    
    
    for i in range(0, len(flat_x)):
        point = Point(flat_x[i], flat_y[i])
        if point.within(poly) == False:
            flat_z[i] = nan
    
    # reshape xyz
    #grid_x = flat_x.reshape(shape(grid_x))
    #grid_y = flat_y.reshape(shape(grid_y))
    grid_z = flat_z.reshape(shape(grid_z))
    
    #grid_x = grid_x[:,::-1]
    #grid_y = grid_y[:,::-1]
    grid_z = grid_z[:,::-1]
    
    '''
    write mesh
    '''
    print 'Writing', key, 'geoTIFF...'
    output_file = ''.join(('geoTIFF', sep, 'GSC2015_hazgrd_', period, '_', key, '.tiff'))
    
    # Create gtif
    nbands = 1
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_file, nlon, nlat, nbands, gdal.GDT_Float32)
    
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    #dst_ds.SetGeoTransform([ minlon, 0.1, 0, minlat, 0, 0.1 ] ) - this works, but map is upsidedown
    dst_ds.SetGeoTransform([ minlon, resolution, 0, maxlat, 0, -resolution])
      
    # set the reference info 
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS("WGS84") # alt NAD27, NAD83
    #srs.SetLCC() # to set Lambert Conformal Conic
    dst_ds.SetProjection(srs.ExportToWkt())
    
    # write the band
    dst_ds.GetRasterBand(1).WriteArray(grid_z.T)
    dst_ds = None # to close file

# testing
'''
src_ds = gdal.Open('NBCC2015_Canada_PGVall_P0.000404.tiff')
src_ds.GetGeoTransform()
srcband = src_ds.GetRasterBand(1)
srcband.GetMetadata()
'''

