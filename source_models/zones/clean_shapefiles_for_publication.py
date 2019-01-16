# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 17:01:44 2018

@author: u56903
"""

import shapefile
from sys import argv
from os import path
from numpy import array, unique, where
from tools.nsha_tools import get_field_data
from mapping_tools import get_WGS84_area
from shapely.geometry import Polygon
from tools.source_shapefile_builder import get_preferred_catalogue, \
                                           get_completeness_model, get_aus_shmax_vectors, \
                                           get_rate_adjust_factor, build_source_shape, \
                                           get_ul_seismo_depths, get_neotectonic_domain_params, \
                                           aggregate_intraslab_sources

shpPath = argv[1]

# read shapefile
print 'Reading source shapefile...'
sf = shapefile.Reader(shpPath)
shapes = sf.shapes()
records = sf.records()
fields = sf.fields

# get polygons
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

'''
newFields = []
fieldType = []
fieldSize = []
fieldDecimal = []
for f in fields[1:]:
    newFields.append(f[0])
    fieldType.append(f[1])
    fieldSize.append(f[2])
    fieldDecimal.append(f[3])
'''
##############################################################################
# build Mmax weights dictionary
##############################################################################

# get Mmax dict from source model weights
ssccsv = '../../shared/seismic_source_model_weights_rounded_p0.4.csv'

lines = open(ssccsv).readlines()

trt = []
mxv = []
mxw = []
for line in lines:
    dat = line.strip().split(',')
    if dat[1] == 'Mmax':
        trt.append(dat[2])
        mxv.append(float(dat[3]))
        mxw.append(float(dat[4]))
        
        # use Proterozoic as generic "Cratonic"
        if dat[2] == 'Proterozoic':
            trt.append('Cratonic')
            mxv.append(float(dat[3]))
            mxw.append(float(dat[4]))
            

trt = array(trt)  
mxv = array(mxv)
mxw = array(mxw)      
unq_trt = unique(trt)

# now make list of dicts
mx_dict = {}
for ut in unq_trt:
    idx = where(trt == ut)[0]
    td = {'trt':ut, 'mx_vals':mxv[idx], 'mx_wts':mxw[idx]}    
    mx_dict[ut] = td
    #print td

##############################################################################
# set shapefile fields
##############################################################################

oldFields = ['SRC_NAME', 'CODE', 'SRC_TYPE', 'CLASS', 'SRC_WEIGHT', 'RTE_ADJ_F',
               'DEP_BEST', 'DEP_UPPER', 'DEP_LOWER', 'USD', 'LSD', 'OW_LSD',
               'MIN_MAG', 'MIN_RMAG', 'MMAX_BEST', 'MMAX_LOWER', 'MMAX_UPPER',
               'N0_BEST', 'N0_LOWER', 'N0_UPPER', 'BVAL_BEST', 'BVAL_LOWER',
               'BVAL_UPPER', 'BVAL_FIX', 'BVAL_FIX_S', 'YCOMP', 'MCOMP',
               'CAT_YMAX', 'PREF_STK', 'PREF_DIP', 'PREF_RKE', 'SHMAX',
               'SHMAX_SIG', 'TRT', 'GMM_TRT', 'DOMAIN', 'CAT_FILE']

newFields = ['SRC_NAME','CODE','DOMAIN','CLASS', 'AREA', 'DEP_BEST', \
             'DEP_UPPER','DEP_LOWER','USD','LSD','OW_LSD','MIN_MAG','MIN_RMAG', \
             'MMAX','MMAX_WGTS','N0_BEST','N0_LOWER','N0_UPPER','RTE_ADJ_F','BVAL_BEST', \
             'BVAL_LOWER','BVAL_UPPER','BVAL_FIX','BVAL_FIX_S','PREF_STK', \
             'PREF_DIP','PREF_RKE','SHMAX','SHMAX_SIG','SCALING_REL', \
             'ASP_RATIO','TRT','GMM_TRT','YCOMP','MCOMP','CAT_YMAX','CAT_FILE']

csvLines = ','.join(newFields) + '\n'

fieldType = ['C','C','F','F','F','F','F','F','F','F','F','F','F','C','C','F', \
             'F','F','F','F','F','F','F','F','F','F','F','F','F','C','F','C','C','C','C','F','C']

fieldSize = [100,12,2,4,10,6,6,6,6,6,5,4,4,30,30,8,8,8,6,6,6,6,6,6,5,5,5,5,5,20,3,30,30,70,50,8,50]

fieldDecimal = [0,0,0,1,0,1,1,1,1,1,1,2,2,0,0,5,5,5,4,3,3,3,3,3,1,1,1,1,1,0,1,0,0,0,0,3,0]

domains = get_field_data(sf, 'DOMAIN', 'float')
trts = get_field_data(sf, 'TRT', 'str')
old_mmaxs = get_field_data(sf, 'MMAX_BEST', 'float')
##############################################################################
# build shapefile
##############################################################################

# set shapefile fields
w = shapefile.Writer(shapefile.POLYGON)
for fn, ft, fs, fd in zip(newFields, fieldType, fieldSize, fieldDecimal):
    if ft == 'F':
        w.field(fn, ft, fs, fd)
    else:
        w.field(fn, ft, str(fs))
    
# now loop through records and get data
i = 0
for record, shape, poly in zip(records, shapes, polygons):
    
    # loop through fields and populate data 
    data = []
    for nf, ft in zip(newFields, fieldType):
        # first use old records
        for of in oldFields:
            if nf == of:
                if ft == 'C':
                    value = get_field_data(sf, nf, 'str')[i]
                else:
                    value = get_field_data(sf, nf, 'float')[i]
                    
        # check if new field and add new data
        domain = domains[i]
        trt = trts[i]
        old_mmax = old_mmaxs[i]
        
        if nf == 'SCALING_REL':            
            if domain <= 7.:
                value = 'Leonard2014_SCR'
            elif domain == 8 or domain == 9:
                value = 'WC1994'
            elif domain == 10:
                value = 'StrasserInterface'
            elif domain == 11:
                value = 'StrasserIntraslab'
                
        elif nf == 'ASP_RATIO':
            if domain <= 7.:
                value = 1.5 # balance between L14 and Cea14 surface rupture lengths
            elif domain == 8 or domain == 9:
                value = 1.5
            elif domain == 10:
                value = 1.5 # based on approx AH interface apect ratios at Mw 8
            elif domain == 11:
                value = 1.2 # based on approx AH intraslab apect ratios at Mw 7.5
                
        elif nf == 'MIN_MAG':            
            if domain <= 7:
                value = 4.5
            elif domain == 8 or domain == 9:
                value = 5.5
            elif domain == 10:
                value = 6.5
            elif domain == 11:
                value = 5.5

        elif nf == 'MMAX':
            try:
                mxArray = mx_dict[trt]['mx_vals']
            except:
                mxArray = [old_mmax-0.2, old_mmax-0.1, old_mmax, old_mmax+0.1, old_mmax+0.2]
                
            mxStr = [str(x) for x in mxArray]
            value  = ';'.join(mxStr)
        
        elif nf == 'MMAX_WGTS':
            try:
                mxwArray  = mx_dict[trt]['mx_wts']
            except:
                mxwArray  = [0.02, 0.14, 0.68, 0.14, 0.02]
            
            mxwStr = [str(x) for x in mxwArray]
            value  = ';'.join(mxwStr)
        
        elif nf == 'OW_LSD':
            value *= -1.
            
        elif nf == 'AREA':
            value = round(get_WGS84_area(poly))
        
        # add data to record
        data.append(value)
    
    print 'Writing record:', i
    
    # only write continental zones
    if domain < 9:
        # set shape polygon
        w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
        # write new records
        w.record(data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7], \
                 data[8],data[9],data[10],data[11],data[12],data[13],data[14], \
                 data[15],data[16],data[17],data[18],data[19],data[20],data[21], \
                 data[22],data[23],data[24],data[25],data[26],data[27],data[28], \
                 data[29],data[30],data[31],data[32],data[33],data[34],data[35],data[36])
    
        # make csv file while I'm at it
        csvLines += ','.join((str(data[0]),str(data[1]),str(data[2]),str(data[3]),str(data[4]),str(data[5]),str(data[6]),str(data[7]), \
                              str(data[8]),str(data[9]),str(data[10]),str(data[11]),str(data[12]),str(data[13]),str(data[14]), \
                              str(data[15]),str(data[16]),str(data[17]),str(data[18]),str(data[19]),str(data[20]),str(data[21]), \
                              str(data[22]),str(data[23]),str(data[24]),str(data[25]),str(data[26]),str(data[27]),str(data[28]), \
                              str(data[29]),str(data[30]),str(data[31]),str(data[32]),str(data[33]),str(data[34]),str(data[35]),str(data[36]))) + '\n'
        
    # increment by 1
    i += 1 
            
# now save area shapefile
outfile = path.split(shpPath)[-1]
outshp = outfile.replace('NSHA18_MFD.shp', 'SSM.shp')
outpath = path.join('shapefiles', 'Publication', outshp)
w.save(outpath)

prjfile = outpath.strip().split('.shp')[0]+'.prj'
f = open(prjfile, 'wb')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()

# now write csv
outcsv = outshp[:-4]+'_Summary.csv'
outpath = path.join('shapefiles', 'Publication', outcsv)
f = open(outpath, 'wb')
f.write(csvLines)
f.close()



























