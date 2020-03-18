7
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:37:07 2018

@author: u56903
"""

def get_completeness_model(src_codes, src_shapes, domains, singleCorner):
    '''
    singleCorner
        1 = do singleCorner (True)
        0 = do not do singleCorner (False)
    '''
    
    from os import path
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    from mapping_tools import distance
    
    # load completeness shp
    if singleCorner == 1:
        compshp = path.join('..','Other','Mcomp_NSHA18_single.shp') # single corner 
    else:
        #compshp = path.join('..','Other','Mcomp_NSHA18_multi.shp') # multi corner 
        compshp = path.join('..','Other','gridded_polygons_3d_completeness.shp') # gridded model for updated Mc - Jan 2020
    
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
        # get centroid of completeness sources
        clon, clat = get_shp_centroid(poly.points)
        point = Point(clon, clat)
        print(clon, clat)
        
        # loop through target and find point in poly    
        mccompFound = False
        dist_to_comp_cent = 9999.
        for i, mc_shape in enumerate(mc_shapes):
            mc_poly = Polygon(mc_shape.points)
            mclon, mclat = get_shp_centroid(mc_shape.points)
            
            # check if target centroid in completeness poly
            if point.within(mc_poly) or point.touches(mc_poly):
                # get dist to centroids
                rngkm = distance(clat, clon, mclat, mclon)[0]
                if rngkm < dist_to_comp_cent:
                    tmp_ycmp = mc_ycomp[i]
                    tmp_mcmp = mc_mcomp[i]
                    mccompFound = True
                    
        # now fill completeness if True
        if mccompFound == True:
            ycomp.append(tmp_ycmp)
            mcomp.append(tmp_mcmp)
        
        # if no Mcomp model assigned, use conservative model
        elif mccompFound == False:
            if dom >= 1 and dom <= 8:
                # for single-corner
                if singleCorner == 1:
                    ycomp.append('1980;1980')
                    mcomp.append('3.5;3.5')
            
                # for mult-corner
                else:
                    ycomp.append('1980;1964;1900')
                    mcomp.append('3.5;5.0;6.0')
                                
            # use approx ISC-GEM completeness
            else:
                ycomp.append('1975;1964;1904')
                mcomp.append('5.75;6.25;7.5')
            
        # set rmin range
        min_rmag.append(max([3.0, float(mcomp[-1].split(';')[0])]))
        
    return ycomp, mcomp, min_rmag

def get_completeness_model_vertex(src_codes, src_shapes, domains, singleCorner):
    '''
    singleCorner
        1 = do singleCorner (True)
        0 = do not do singleCorner (False)
    '''
    
    from os import path
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    
    # load completeness shp
    if singleCorner == 1:
        compshp = path.join('..','Other','Mcomp_NSHA18_single.shp') # single corner 
    else:
        #compshp = path.join('..','Other','Mcomp_NSHA18_multi.shp') # multi corner 
        compshp = path.join('..','Other','gridded_polygons_3d_completeness.shp') # gridded model for updated Mc - Jan 2020
    
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
        tmp_mcomp = '-99;-99' # dummy value
        mccompFound = False
                
        # get centroid of completeness sources
        for clon, clat in poly.points:
            point = Point(clon, clat)
            
            # loop through target and find point in poly    
            for i in range(0, len(mc_shapes)):
                mc_poly = Polygon(mc_shapes[i].points)
                
                # check if target vertex in completeness poly
                if point.within(mc_poly): 
                    # select most conservative option
                    if float(mc_mcomp[i].strip().split(';')[0]) > float(tmp_mcomp.strip().split(';')[0]):
                        tmp_ycomp = mc_ycomp[i]
                        tmp_mcomp = mc_mcomp[i]
                        mccompFound = True
            
        # if no Mcomp model assigned, use conservative model
        if mccompFound == False:
            if dom <= 8:
                # for single-corner
                if singleCorner == 1:
                    tmp_ycomp = '1980;1980'
                    tmp_mcomp = '3.5;3.5'
            
                # for mult-corner
                else:
                    tmp_ycomp = '1980;1964;1900'
                    tmp_mcomp = '3.5;5.0;6.0'
                                
            # use approx ISC-GEM completeness
            else:
                tmp_ycomp = '1975;1964;1904'
                tmp_mcomp = '5.75;6.25;7.5'
                
            
        ycomp.append(tmp_ycomp)
        mcomp.append(tmp_mcomp)
        
        # set rmin range
        min_rmag.append(max([3.0, float(mcomp[-1].split(';')[0])]))
        
    return ycomp, mcomp, min_rmag

def get_completeness_model_point(clat, clon, singleCorner):
    '''
    singleCorner
        1 = do singleCorner (True)
        0 = do not do singleCorner (False)
        
        assume AU, dom = 0
    '''
    dom = 0
    from os import path, getcwd
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    from mapping_tools import distance
    
    # load completeness shp
    if getcwd().startswith('/Users'):
        if singleCorner == 1:
            compshp = path.join('/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Other/Mcomp_NSHA18_single.shp') # single corner 
        else:
            compshp = path.join('/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Other/gridded_polygons_3d_completeness.shp') # multi corner 
    else:
        if singleCorner == 1:
            compshp = path.join('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/trev/NSHA2018/source_models/zones/shapefiles/Other/Mcomp_NSHA18_single.shp') # single corner 
        else:
            #compshp = path.join('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/trev/NSHA2018/source_models/zones/shapefiles/Other/Mcomp_NSHA18_multi.shp') # multi corner 
            compshp = path.join('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/trev/NSHA2018/source_models/zones/shapefiles/Other/Mcomp_NSHA18_multi_20191217.shp') # multi corner 
    
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
    point = Point(clon, clat)
    print(clon, clat)
        
    # loop through target and find point in poly    
    mccompFound = False
    dist_to_comp_cent = 9999.
    for i, mc_shape in enumerate(mc_shapes):
        mc_poly = Polygon(mc_shape.points)
        mclon, mclat = get_shp_centroid(mc_shape.points)
        
        # check if target centroid in completeness poly
        if point.within(mc_poly) or point.touches(mc_poly):
            # get dist to centroids
            rngkm = distance(clat, clon, mclat, mclon)[0]
            print(rngkm)
            if rngkm < dist_to_comp_cent:
                ycomp.append(mc_ycomp[i])
                mcomp.append(mc_mcomp[i])
                mccompFound = True
                print(mc_poly)
        
    # if no Mcomp model assigned, use conservative model
    if mccompFound == False:
        if dom <= 8:
            # for single-corner
            if singleCorner == 1:
                ycomp = '1980;1980'
                mcomp = '3.5;3.5'
        
            # for mult-corner
            else:
                ycomp = '1980;1964;1900'
                mcomp = '3.5;5.0;6.0'
                            
        # use approx ISC-GEM completeness
        else:
            ycomp = '1975;1964;1904'
            mcomp = '5.75;6.25;7.5'
        
    # set rmin range
    min_rmag.append(max([3.0, float(mcomp.split(';')[0])]))
        
    return ycomp, mcomp, min_rmag
    
# need to ensure upper/lower seismo depths consistent with Domains edits
def get_ul_seismo_depths(target_codes, target_usd, target_lsd):
    from os import path
    import shapefile
    from numpy import array, median, std
    from tools.nsha_tools import get_field_data
    
    shmaxshp = path.join('..','Domains','Domains_NSHA18_single_Mc.shp')

    print('Reading SHmax shapefile...')
    sf = shapefile.Reader(shmaxshp)
        
    # get shmax attributes
    source_codes = get_field_data(sf, 'CODE', 'str')
    source_usd = get_field_data(sf, 'USD', 'float')
    source_lsd = get_field_data(sf, 'LSD', 'float')
    
    for j, tc in enumerate(target_codes):
        matchCodes = False
        for i, sc in enumerate(source_codes):
            if tc == sc:
                target_usd[j] = source_usd[i]
                target_lsd[j] = source_lsd[i]
                matchCodes = True
        
        # no match
        if matchCodes == False:
           print('  Cannot match seis depths:', tc)
           
    return target_usd, target_lsd
    
    
# get neotectonic domain number and Mmax from zone centroid
def get_neotectonic_domain_params(target_sf, target_trt, refShpFile):
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    
    # load target shapefile
    polygons = target_sf.shapes()
    
    # load domains shp
    domshp = open(refShpFile).read()
    dsf = shapefile.Reader(domshp)
    
    # get domains
    neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
    neo_min_reg = get_field_data(dsf, 'MIN_RMAG', 'float')
    neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
    neo_bval = get_field_data(dsf, 'BVAL_BEST', 'float')
    neo_bval_l = get_field_data(dsf, 'BVAL_LOWER', 'float')
    neo_trt  = get_field_data(dsf, 'TRT', 'str')
    
    # get bval sigma
    bval_sig = neo_bval_l - neo_bval
    
    # get domain polygons
    dom_shapes = dsf.shapes()
    domain = []
    min_rmag = []
    mmax = []
    trt = []
    bval_fix = []
    bval_sig_fix = []
    
    # loop through target zones
    for poly, ttrt in zip(polygons, target_trt):
        # get centroid of target sources
        clon, clat = get_shp_centroid(poly.points)
        point = Point(clon, clat)
        
        # loop through domains and find point in poly
        matchidx = -99
        for i in range(0, len(dom_shapes)):
            # make sure trts match
            if neo_trt[i] == ttrt:
                dom_poly = Polygon(dom_shapes[i].points)
                
                # check if target centroid in domains poly
                if point.within(dom_poly):
                    matchidx = i
                
        # set dummy values
        if matchidx == -99:
            domain.append(-99)
            min_rmag.append(3.5)
            mmax.append(-99)
            trt.append('')
            bval_fix.append(-99)
            bval_sig_fix.append(-99)
        # fill real values
        else:
            domain.append(neo_doms[matchidx])
            min_rmag.append(neo_min_reg[matchidx])
            mmax.append(neo_mmax[matchidx])
            trt.append(neo_trt[matchidx])
            bval_fix.append(neo_bval[matchidx])
            bval_sig_fix.append(bval_sig[matchidx])
        
    return domain, min_rmag, mmax, trt, bval_fix, bval_sig_fix

# get neotectonic domain number and Mmax from zone centroid
def get_simple_neotectonic_domain_params(target_sf, refShpFile):
    import shapefile
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data, get_shp_centroid
    
    # load target shapefile
    polygons = target_sf.shapes()
    
    # load domains shp
    domshp = open(refShpFile).read()
    dsf = shapefile.Reader(domshp)
    
    # get domains
    neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
    neo_min_reg = get_field_data(dsf, 'MIN_RMAG', 'float')
    neo_mmax = get_field_data(dsf, 'MMAX_BEST', 'float')
    neo_bval = get_field_data(dsf, 'BVAL_BEST', 'float')
    neo_bval_l = get_field_data(dsf, 'BVAL_LOWER', 'float')
    neo_trt  = get_field_data(dsf, 'TRT', 'str')
    neo_usd  = get_field_data(dsf, 'USD', 'float')
    neo_lsd  = get_field_data(dsf, 'LSD', 'float')
    neo_dep_b  = get_field_data(dsf, 'DEP_BEST', 'float')
    neo_dep_u  = get_field_data(dsf, 'DEP_UPPER', 'float')
    neo_dep_l  = get_field_data(dsf, 'DEP_LOWER', 'float')
    
    # get bval sigma
    bval_sig = neo_bval_l - neo_bval
    
    # get domain polygons
    dom_shapes = dsf.shapes()
    domain = []
    min_rmag = []
    mmax = []
    trt = []
    bval_fix = []
    bval_sig_fix = []
    usd = []
    lsd = []
    dep_b = []
    dep_u = []
    dep_l = []
    
    # loop through target zones
    for poly in polygons:
        # get centroid of target sources
        clon, clat = get_shp_centroid(poly.points)
        point = Point(clon, clat)
        
        # loop through domains and find point in poly
        matchidx = -99
        for i in range(0, len(dom_shapes)):
            # make sure trts match
            dom_poly = Polygon(dom_shapes[i].points)
                
            # check if target centroid in domains poly
            if point.within(dom_poly):
                matchidx = i
                
        # set dummy values
        if matchidx == -99:
            domain.append(-99)
            min_rmag.append(3.5)
            mmax.append(8.5)
            trt.append('')
            bval_fix.append(-99)
            bval_sig_fix.append(-99)
            usd.append(-99)
            lsd.append(-99)
            dep_b.append(-99)
            dep_u.append(-99)
            dep_l.append(-99)
        # fill real values
        else:
            domain.append(neo_doms[matchidx])
            min_rmag.append(neo_min_reg[matchidx])
            mmax.append(neo_mmax[matchidx])
            trt.append(neo_trt[matchidx])
            bval_fix.append(neo_bval[matchidx])
            bval_sig_fix.append(bval_sig[matchidx])
            usd.append(neo_usd[matchidx])
            lsd.append(neo_lsd[matchidx])
            dep_b.append(neo_dep_b[matchidx])
            dep_u.append(neo_dep_u[matchidx])
            dep_l.append(neo_dep_l[matchidx])
    
    return domain, min_rmag, mmax, trt, bval_fix, bval_sig_fix, usd, lsd, dep_b, dep_u, dep_l

# use Rajabi_2016 shmax vectors - gets median & std within a source zone    
def get_aus_shmax_vectors(src_codes, src_shapes):
    from os import path
    import shapefile
    from numpy import array, median, std
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data
    
    print('Reading SHmax shapefile...')
    try:
        #shmaxshp = path.join('..','Other','SHMax_Rajabi_2016.shp')
        shmaxshp = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/trev/NSHA2018/source_models/zones/shapefiles/Other/SHMax_Rajabi_2016.shp'
        sf = shapefile.Reader(shmaxshp)
    except:
        shmaxshp = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Other/SHMax_Rajabi_2016.shp'
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
            
            # check sigma and make sure it is at least +/- 15 degrees
            shmax_sig.append(max([std(array(shm_in)), 15.]))        
            
            print('Getting SHmax for', code)
        
        # if no points in polygons, get nearest neighbour
        else:
            print('Getting nearest neighbour...')
            min_dist = 9999.
            for shmlo, shmla, shm in zip(shmax_lon, shmax_lat, shmax):
                pt = Point(shmlo, shmla)
                pt_dist = pt.distance(Polygon(poly.points))
                if pt_dist < min_dist:
                    min_dist = pt_dist
                    shm_near = shm
            
            shmax_pref.append(shm_near) # set nearest neighbour
            shmax_sig.append(15.) # set std manually
            
    return shmax_pref, shmax_sig
    
# use gridded b-value to get fixed b
def get_gridded_bvalue(src_codes, src_shapes):
    from os import path
    from numpy import array, nanmedian, std, loadtxt, isnan
    from shapely.geometry import Point, Polygon
    from tools.nsha_tools import get_field_data
    
    print('Reading b-value grid...')
    try:
        #shmaxshp = path.join('..','Other','SHMax_Rajabi_2016.shp')
        bgrid_csv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/trev/NSHA2018/source_models/zones/shapefiles/Other/SHMax_Rajabi_2016.shp'
        # get shmax attributes
        blons, blats, bvals, bsigs = loadtxt(bgrid_csv, delimiter=',', skiprows=1)
    
    except:
        bgrid_csv = '/Users/trev/Documents/Manuscripts/manuscripts/conferences/2020_wcee/gridded_bval.csv'
        # get shmax attributes
        data = loadtxt(bgrid_csv, delimiter=',', skiprows=1)
    
    blons = data[:, 0]
    blats = data[:, 1]
    bvals = data[:, 2]
    bsigs = data[:, 3]
    
    ###############################################################################
    # get preferred strike
    ###############################################################################
    bval_fix = []
    bval_sig_fix  = []
    
    for code, poly in zip(src_codes, src_shapes):
        # get shmax points in polygon
        bval_in = []
        bval_sig_in = []
        
        # now loop through earthquakes in cat
        for blo, bla, bv, bs in zip(blons, blats, bvals, bsigs):
            
            # check if pt in poly and compile mag and years
            pt = Point(blo, bla)
            if pt.within(Polygon(poly.points)):
                bval_in.append(bv)
                bval_sig_in.append(bs)
        
        if len(bval_in) > 0: 
            bval_fix.append(nanmedian(array(bval_in)))
            
            # check sigma and make sure it is at least +/- 15 degrees
            bval_sig_fix.append(nanmedian(array(bval_sig_in)))
            
            print('Getting b-value for', code)
        
        # if no points in polygons, get nearest neighbour
        else:
            print('Getting nearest neighbour...')
            min_dist = 9999.
            for blo, bla, bv, bs in zip(blons, blats, bvals, bsigs):
                pt = Point(blo, bla)
                pt_dist = pt.distance(Polygon(poly.points))
                if pt_dist < min_dist:
                    min_dist = pt_dist
                    bv_near = bv
                    bs_near = bs
            
            if isnan(bv_near) or bv_near < 0:
                bv_near = 1.0
                bs_near = 0.15
                
            bval_fix.append(bv_near) # set nearest neighbour
            bval_sig_fix.append(bs_near) # set std manually
            
    return bval_fix, bval_sig_fix 
    

def get_rate_adjust_factor(newshp, newField, origshp, origField):
    import shapefile
    from numpy import ones
    from tools.mapping_tools import get_WGS84_area
    from shapely.geometry import Polygon
    from tools.nsha_tools import get_field_data
    
    print('\nChecking shape geometrties...')
    
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
                    
                    print('    ',newCode,'rate adjustment factor:',rte_adj_fact[i])
                    
        i += 1
        
    return rte_adj_fact
    
# gets preferred catalogue - uses NSHA cat if all vertices inside GG_cat_polygon
def get_preferred_catalogue(targetshpfile):
    import shapefile
    from shapely.geometry import Point, Polygon
    #from tools.nsha_tools import get_field_data
    
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
                tmpcat = 'ISC-GEM_V5_hmtk_GK74_declustered_clip.csv'
                
        # now append catalogue
        cat.append(tmpcat)
        
    return cat
    
# function to aggregate different intraslab sources to obtian improved earthquake recurrence
def aggregate_intraslab_sources(src_code, src_class):
    from os import path
    
    # set aggregation code file relative to shapfiles
    aggFile = path.join('..','..','..','banda','banda_deep_sources_NSHA18_aggregation_map.csv')
    
    # read file
    lines = open(aggFile).readlines()[1:]
    
    # loop thru lines and find matching code
    for line in lines:
        dat = line.strip().split(',')
        agg_code = dat[2]
        agg_class = float(dat[-1])
        
        # if match, reset zone class
        if agg_code == src_code:
            src_class = agg_class
    
    return src_class

def build_source_shape(outshp, src_shapes, src_names, src_codes, zone_class, \
                       rte_adj_fact, dep_b, dep_u, dep_l, usd, lsd, \
                       min_rmag, mmax, bval_fix, bval_sig_fix, \
                       ycomp, mcomp, pref_stk, pref_dip, pref_rke, \
                       shmax_pref, shmax_sig, trt, dom, prefCat):
                           
    import shapefile
    from numpy import array, ones_like, where
    
    # many eqs within Aust get left out if LSD too conservative    
    overwrite_lsd = 999 * ones_like(lsd)
    
    #idx = array(dom) <= 8 # hardwire for continental sources
    #overwrite_lsd[idx] = 999 # in km
    
    # set overwright_lsd for insalb sources
    idx = where(array(dom) == 11)[0]
    overwrite_lsd[idx] = lsd[idx] + 200 # in km
    
    # get effective trt for GMM assignment
    gmm_trt = []
    for t in trt:
        if t == 'Cratonic':
            gmm_trt.append(t)
        elif t == 'Active' or t == 'Extended' or t == 'Non_cratonic' or t == 'Oceanic':
            gmm_trt.append('Non_cratonic')
        elif t == 'Intraslab' or t == 'Interface':
            gmm_trt.append('Subduction')
        else:
            gmm_trt.append('Non_cratonic')
    
    # set shapefile to write to
    try:
        w = shapefile.Writer(shapefile.POLYGON)
    except:
        w = shapefile.Writer(outshp[:-4], shapeType=5)
    w.field('SRC_NAME','C','100')
    w.field('CODE','C','12')
    w.field('SRC_TYPE','C','10')
    w.field('CLASS','C','10')
    w.field('SRC_WEIGHT','F', 8, 2)
    w.field('RTE_ADJ_F','F', 6, 4)
    w.field('DEP_BEST','F', 6, 1)
    w.field('DEP_UPPER','F', 6, 1)
    w.field('DEP_LOWER','F', 6, 1)
    w.field('USD','F', 4, 1)
    w.field('LSD','F', 4, 1)
    w.field('OW_LSD','F', 4, 1)
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
    w.field('CAT_YMAX', 'F', 8, 3)
    w.field('PREF_STK','F', 6, 2)
    w.field('PREF_DIP','F', 6, 2)
    w.field('PREF_RKE','F', 6, 2)
    w.field('SHMAX','F', 6, 2)
    w.field('SHMAX_SIG','F', 6, 2)
    w.field('TRT','C','100')
    w.field('GMM_TRT','C','100')
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
    cat_ymax = -99
        
    # loop through original records
    for i, shape in enumerate(src_shapes):
    
        # set shape polygon
        try:
            w.line(parts=[shape.points], shapeType=shapefile.POLYGON)
        except:
            w.poly([shape.points])
            
        # write new records
        if i >= 0:
            #print(src_names[i], src_codes[i], src_ty, zone_class[i], src_wt, rte_adj_fact[i], dep_b[i], dep_u[i], dep_l[i], usd[i], lsd[i])
            #print(overwrite_lsd[i], min_mag, min_rmag[i], mmax[i], mmax[i]-0.2, mmax[i]+0.2, n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix[i], bval_sig_fix[i])
            #print(
            w.record(src_names[i], src_codes[i], src_ty, zone_class[i], src_wt, rte_adj_fact[i], dep_b[i], dep_u[i], dep_l[i], usd[i], lsd[i], \
                     overwrite_lsd[i], min_mag, min_rmag[i], mmax[i], mmax[i]-0.2, mmax[i]+0.2, n0, n0_l, n0_u, bval, bval_l, bval_u, bval_fix[i], bval_sig_fix[i], \
                     ycomp[i], mcomp[i], cat_ymax, pref_stk[i], pref_dip[i], pref_rke[i], shmax_pref[i], shmax_sig[i], trt[i], gmm_trt[i], dom[i], prefCat[i])
            
    # now save area shapefile
    try:
        w.save(outshp)
    except:
        w.close()
    
    # write projection file
    print(outshp)
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'w')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()