# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:07:46 2015

@author: tallen
"""
# https://gis.stackexchange.com/questions/130199/changing-color-of-raster-images-based-on-their-data-values-gdal
#from pyrate.shared import Ifg, DEM
#import subprocess
import sys
import os
import tempfile
import numpy as np

def gen_color_file(input_file):
    fp, temp_file = tempfile.mkstemp(suffix='.txt')

    dem = DEM(input_file)
    dem.open()
    phase_data = dem.height_band.ReadAsArray()

    max_ph = np.nanmax(phase_data)
    min_ph = np.nanmin(phase_data)
    range_ph = max_ph-min_ph
    colors = ['black', 'blue', 'yellow', 'orange', 'red', 'white']
    with open(temp_file, 'w') as f:
        for i, c in enumerate(colors[:-1]):
            f.write(str(int(min_ph + (i + 1)*range_ph/len(colors))) +
                    ' ' + c + '\n')
        f.write(str(int(max_ph - range_ph/len(colors))) +
                ' ' + colors[-1] + '\n')
    os.close(fp)
    return temp_file
'''
cmd = "gdaldem color-relief " + input_file \
      + ' ' + color_file + ' ' + output_file
      subprocess.check_call(cmd, shell=True)
'''    
'''
# read sol file
'''
from sys import argv
from os import sep, path, mkdir, system, getcwd
from numpy import array, mgrid, nan, shape, hstack, isinf, log, exp, interp
from scipy.interpolate import griddata
#from matplotlib.mlab import griddata
from osgeo import osr, gdal
import shapefile	
from shapely.geometry import Point, Polygon
from tools.oq_tools import return_annualised_haz_curves
#from misc_tools import dictlist2array

maxlat = -90
minlat = 90
maxlon = -180
minlon = 180

##############################################################################
# parse hazard grid
##############################################################################

hazCurveGridFile = argv[1]

# check to see if geotiff folder exists
if path.isdir('geotiff') == False:
    mkdir('geotiff')
    
# make out geoTIFF    
period = hazCurveGridFile.split(sep)[-2].split('_')[-1]

# parse grid file
gridDict, imls, investigation_time = return_annualised_haz_curves(hazCurveGridFile)

##############################################################################
# interpolate hazard curve and fill dictionary
##############################################################################

# set grid return periods for hazard curve
probs = array([0.02,0.01375,0.01,0.00445,0.002,0.0021,0.001,0.0005,0.000404,0.0002,0.0001])

grddict = []
alon = []
alat = []

for site in gridDict:
    interpHaz = exp(interp(log(probs[::-1]), log(site['poe_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    
    # fill a temp dictionary
    tmpdict = {'lon':site['lon'], 'lat':site['lat']}
    
    for p, ih in zip(probs, interpHaz):
        if ih < 1E-20:
            tmpdict['P'+str(p)] = 1E-20
        else:
            tmpdict['P'+str(p)] = ih
    
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

##############################################################################
# make mesh
##############################################################################

resolution = 0.025 # degrees
invres = int(1./resolution)
keys = ['P0.0021', 'P0.000404'] # probabilities
pc50 = ['0.1', '0.02']
for key, p50 in zip(keys, pc50):
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
    if getcwd().startswith('/nas'):
        inshape = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/shapefiles//au_maritime_boundary_digitised.shp'
    
    sf = shapefile.Reader(inshape)
    sf = sf.shapes()
    poly = Polygon(sf[0].points)
    
    flat_x = grid_x.flatten()
    flat_y = grid_y.flatten()
    flat_z = grid_z.flatten()
    
    # now loop through points
    print 'Cropping points...'
    for i in range(0, len(flat_x)):
        point = Point(flat_x[i], flat_y[i])
        if point.within(poly) == False:
            flat_z[i] = nan
    
    # reshape xyz
    grid_z = flat_z.reshape(shape(grid_z))
    grid_z = grid_z[:,::-1]
    
    '''
    write mesh
    '''
    print 'Writing', key, 'geoTIFF...'
    output_file = path.join('geotiff', '_'.join(('nsha18',period, p50+'.tiff')))
    
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
    
    '''
    # set color
    ct = gdal.ColorTable()
    # Some examples
    ct.SetColorEntry( 0, (0, 0, 0, 255) )
    ct.SetColorEntry( 1, (0, 255, 0, 255) )
    ct.SetColorEntry( 2, (255, 0, 0, 255) )
    ct.SetColorEntry( 3, (255, 0, 255, 255) )
    # Set the color table for your band
    dst_ds.GetRasterBand( 1 ).SetRasterColorTable( ct )
    '''
    # write the band
    dst_ds.GetRasterBand(1).WriteArray(grid_z.T)
    dst_ds = None # to close file

# testing

src_ds = gdal.Open(path.join('geotiff', '_'.join(('nsha18',period, p50+'.tiff'))))
'''
src_ds.GetGeoTransform()
srcband = src_ds.GetRasterBand(1)
src_ds.GetMetadata()
'''

'''
band = src_ds.GetRasterBand(1)
ct   = band.GetRasterColorTable()
f    = open("rgb_color.txt", 'w+')    
for i in range(ct.GetCount()):
    sEntry = ct.GetColorEntry(i)
    f.write( "  %3d: %d,%d,%d\n" % ( \
      i, \
      sEntry[0],\
      sEntry[1],\
      sEntry[2]))
'''

'''
format of relief file from: http://blog.mastermaps.com/2012/06/creating-color-relief-and-slope-shading.html

0 110 220 110
900 240 250 160 
1300 230 220 170 
1900 220 220 220
2500 250 250 250 

gdaldem color-relief jotunheimen.tif color_relief.txt jotunheimen_colour_relief.tif
'''

#system('gdaldem color-relief -of VRT input.tif rgb_color.txt rgb_output.vrt')
