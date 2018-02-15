# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:37:07 2018

@author: u56903
"""

def get_completeness_model(src_codes, src_shapes, domains):
    from os import path
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    
    # load domains shp
    compshp = path.join('..','Other','Mcomp_NSHA18_smoothed.shp')
    mcsf = shapefile.Reader(compshp)
    
    # get completeness data
    mc_ycomp = get_field_data(mcsf, 'YCOMP', 'str')
    mc_mcomp = get_field_data(mcsf, 'MCOMP', 'str')
    
    # get completeness polygons
    mc_shapes = mcsf.shapes()
    
    # set empty completeness values
    ycomp = []
    mcomp = []
    min_rmag = []
    
    # loop through Mcomp zones
    for code, poly, dom in zip(src_codes, src_shapes, domains):
        # get centroid of leonard sources
        clon, clat = get_shp_centroid(poly.points)
        point = Point(clon, clat)
        print clon, clat
        
        # loop through domains and find point in poly    
        mccompFound = False
        for i in range(0, len(mc_shapes)):
            dom_poly = Polygon(mc_shapes[i].points)
            
            # check if Domains centroid in completeness poly
            if point.within(dom_poly): 
                ycomp.append(mc_ycomp[i])
                mcomp.append(mc_mcomp[i])
                mccompFound = True
        
        # if no Mcomp model assigned, use conservative model
        if mccompFound == False:
            if dom <= 8:
                ycomp.append('1980;1964;1900')
                mcomp.append('3.5;5.0;6.0')
            
            # use ISC-GEM completeness
            else:
                ycomp.append('1980;1920;1900')
                mcomp.append('5.75;6.25;7.5')
            
        # set rmin range
        min_rmag.append(max([3.0, float(mcomp[-1].split(';')[0])]))
        
    return ycomp, mcomp, min_rmag
    
def get_aus_shmax_vectors(src_codes, src_shapes):
    from os import path
    import shapefile
    from numpy import array, median, std
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data
    
    shmaxshp = path.join('..','Other','SHMax_Rajabi_2016.shp')

    print 'Reading SHmax shapefile...'
    sf = shapefile.Reader(shmaxshp)
        
    # get shmax attributes
    shmax_lat = get_field_data(sf, 'LAT', 'float')
    shmax_lon = get_field_data(sf, 'LON', 'float')
    shmax     = get_field_data(sf, 'SHMAX', 'float')
    
    ###############################################################################
    # get preferred strike
    ###############################################################################
    shmax_pref = []
    shmax_sig  = []
    
    for code, poly in zip(src_codes, src_shapes):
        # get shmax points in polygon
        shm_in = []
        
        # now loop through earthquakes in cat
        for shmlo, shmla, shm in zip(shmax_lon, shmax_lat, shmax):
            
            # check if pt in poly and compile mag and years
            pt = Point(shmlo, shmla)
            if pt.within(Polygon(poly.points)):
                shm_in.append(shm)
        
        if len(shm_in) > 0: 
            shmax_pref.append(median(array(shm_in)))
            shmax_sig.append(std(array(shm_in)))
            print 'Getting SHmax for', code
        
        # if no points in polygons, get nearest neighbour
        else:
            print 'Getting nearest neighbour...'
            min_dist = 9999.
            for shmlo, shmla, shm in zip(shmax_lon, shmax_lat, shmax):
                pt = Point(shmlo, shmla)
                pt_dist = pt.distance(Polygon(poly.points))
                if pt_dist < min_dist:
                    min_dist = pt_dist
                    shm_near = shm
            
            shmax_pref.append(shm_near) # set nearest neighbour
            shmax_sig.append(10.) # set std manually
            
    return shmax_pref, shmax_sig


def get_rate_adjust_factor(newshp, newField, origshp, origField):
    import shapefile
    from numpy import ones
    from tools.mapping_tools import get_WGS84_area
    from shapely.geometry import Polygon
    from tools.nsha_tools import get_field_data
    
    print '\nChecking shape geometrties...'
    
    newsf = shapefile.Reader(newshp)
    new_shapes = newsf.shapes()
    new_codes = get_field_data(newsf, newField, 'str')
    
    origsf = shapefile.Reader(origshp)
    orig_shapes = origsf.shapes()
    orig_codes = get_field_data(origsf, origField, 'str')
    
    # set initial values
    rte_adj_fact = ones(len(new_codes))
    i = 0
    for newCode, newPoly in zip(new_codes, new_shapes):
        newPolyArea = get_WGS84_area(Polygon(newPoly.points))
        
        for origCode, origPoly in zip(orig_codes, orig_shapes):
            if origCode == newCode:
                origPolyArea = get_WGS84_area(Polygon(origPoly.points))
                
                if round(newPolyArea, 4) != round(origPolyArea, 4):
                    rte_adj_fact[i] = round(newPolyArea / origPolyArea, 4)
                    
                    print '    ',newCode,'rate adjustment factor:',rte_adj_fact[i]
                    
        i += 1
        
    return rte_adj_fact
    
# gets preferred catalogue - uses NSHA cat if all vertices inside GG_cat_polygon
def get_preferred_catalogue(targetshpfile):
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data
    
    ###############################################################################
    # parse shapefile
    ###############################################################################
    
    catshpfile = '..//Other//ggcat_extent.shp'
    
    csf = shapefile.Reader(catshpfile)
    catshape = csf.shapes()[0] # only one shape
        
    tsf = shapefile.Reader(targetshpfile)
    targetshapes = tsf.shapes()
    
    ###############################################################################
    # find poly points within other polygons
    ###############################################################################
    
    # loop through target zones
    cat = []
    for poly in targetshapes:
        
        # set default catalogue
        tmpcat = 'NSHA18CAT_V0.1_hmtk_declustered.csv'
        
        # now loop through all points in target shape
        for tlon, tlat in poly.points:
            point = Point(tlon, tlat)
            
            # check if point in catshape
            if point.within(Polygon(catshape.points)) == False:
                tmpcat = 'ISC-GEM_hmtk_declustered.csv'
                
        # now append catalogue
        cat.append(tmpcat)
        
    return cat

def build_source_shape(outshp, src_shapes, src_names, src_codes, zone_class, \
                       rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                       min_rmag, mmax, bval_fix, bval_sig_fix, \
                       ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                       shmax_pref, shmax_sig, trt, dom, prefCat):
                           
    import shapefile
    
    # set shapefile to write to
    w = shapefile.Writer(shapefile.POLYGON)
    w.field('SRC_NAME','C','100')
    w.field('CODE','C','10')
    w.field('SRC_TYPE','C','10')
    w.field('CLASS','C','10')
    w.field('SRC_WEIGHT','F', 8, 2)
    w.field('RTE_ADJ_F','F', 6, 4)
    w.field('DEP_BEST','F', 6, 1)
    w.field('DEP_UPPER','F', 6, 1)
    w.field('DEP_LOWER','F', 6, 1)
    w.field('USD','F', 6, 1)
    w.field('LSD','F', 6, 1)
    w.field('MIN_MAG','F', 4, 2)
    w.field('MIN_RMAG','F', 4, 2)
    w.field('MMAX_BEST','F', 4, 2)
    w.field('MMAX_LOWER','F', 4, 2)
    w.field('MMAX_UPPER','F', 4, 2)
    w.field('N0_BEST','F', 8, 5)
    w.field('N0_LOWER','F', 8, 5)
    w.field('N0_UPPER','F', 8, 5)
    w.field('BVAL_BEST','F', 6, 3)
    w.field('BVAL_LOWER','F', 6, 3)
    w.field('BVAL_UPPER','F', 6, 3)
    w.field('BVAL_FIX','F', 6, 3)
    w.field('BVAL_FIX_S','F', 6, 3)
    w.field('YCOMP','C','70')
    w.field('MCOMP','C','50')
    w.field('PREF_STK','F', 6, 2)
    w.field('PREF_DIP','F', 6, 2)
    w.field('PREF_RKE','F', 6, 2)
    w.field('SHMAX','F', 6, 2)
    w.field('SHMAX_SIG','F', 6, 2)
    w.field('TRT','C','100')
    w.field('DOMAIN','F', 2, 0)
    w.field('CAT_FILE','C','50')
    
    src_wt = 1.0
    src_ty = 'area'
    
    min_mag = 4.5
    n0 = -99
    n0_l = -99
    n0_u = -99
    bval = -99
    bval_l = -99
    bval_u = -99
    
    # loop through original records
    for i, shape in enumerate(src_shapes):
    
        # set shape polygon
        w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
            
        # write new records
        if i >= 0:
            w.record(src_names[i], src_codes[i], src_ty, zone_class[i], src_wt, rte_adj_fact[i], dep_b[i], dep_u[i], dep_l[i], usd[i], lsd[i], \
                     min_mag, min_rmag[i], mmax[i], mmax[i]-0.2, mmax[i]+0.2, n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix, bval_sig_fix, \
                     ycomp[i], mcomp[i], pref_stk[i], pref_dip[i], pref_rke[i], shmax_pref[i], shmax_sig[i], trt[i], dom[i], prefCat[i])
            
    # now save area shapefile
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'wb')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()