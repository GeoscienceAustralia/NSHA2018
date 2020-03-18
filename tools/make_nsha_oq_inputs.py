#####################################################################################
# imports source model from get_NSHA18_mfds.py and outputs an OpenQuake source file with
# collapsed earthquake rates

# Usage: python make_openquake_source_file.py <path to pkl file> <output base name>
#            e.g. python make_openquake_source_file.py SWCan_T3EclC1.pkl swcan
#####################################################################################

def make_collapse_occurrence_text(m, min_mag, binwid, meta, mx_dict):
    from numpy import array, zeros, argmax, zeros_like
    from tools.oq_tools import get_oq_incrementalMFD
    '''
    best_beta: boolean; True if only using best curve
    mx_vals: array of Mmax values
    mx_wts: array of Mmax weights
    '''
    
    # set mx arrays from mx dictionary using zone trt
    # first try EE mags
    try:
        mx_vals = mx_dict[m['trt']]['mx_vals']
        mx_wts  = mx_dict[m['trt']]['mx_wts']
    
    # else use Mmax from model shapefile and theoretical weights
    except:
        mx = m['max_mag'][0]
        mx_vals = [mx-0.2, mx-0.1, mx, mx+0.1, mx+0.2]
        mx_wts  = [0.02, 0.14, 0.68, 0.14, 0.02]
    
    # REMOVE ME WHEN FINISHED TESTING
    #mx_wts = [0.1, 0.2, 0.4, 0.2, 0.1]
        
    # if using one Mx only
    if meta['one_mx'] == True:
        # find highest weighted
        if meta['mx_idx'] == -1:
            # find and set model Mmax
            mxidx = argmax(array(mx_wts))
            
            # reset weight
            mx_wts = zeros_like(mx_wts)
            mx_wts[mxidx] = 1.0
        
        # else use pre-defined index
        else:
            # reset weight
            mx_wts = zeros_like(mx_wts)
            mx_wts[meta['mx_idx']] = 1.0
        
    wtd_list  = []
    maglen = 0
    
    for beta_val, beta_wt, N0 in zip(m['src_beta'], meta['beta_wts'], m['src_N0']):
        # get effective rate
        effN0 = N0 * m['src_weight'] * m['rate_adj_fact']
        
        for mx, mxwt in zip(mx_vals, mx_wts):
            
            betacurve, mrange = get_oq_incrementalMFD(beta_val, effN0, \
                                                      min_mag, mx, \
                                                      binwid)
            
            wtd_list.append(betacurve * beta_wt * mxwt)
            
            # get max length of arrays
            if len(betacurve) > maglen:
                maglen = len(betacurve)

    # sum rates for each branch
    wtd_rates = zeros(maglen)
    for rates in wtd_list:
        # go element by element
        for r, rate in enumerate(rates):
            wtd_rates[r] += rate
    
    #print(wtd_rates[0]    
       
    # convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(wtd_rates[0:-1])):
        occ_rates.append(wtd_rates[b] - wtd_rates[b+1])
    occ_rates.append(wtd_rates[-1])
    
    # make text object                        
    octxt = str('%0.5e' % occ_rates[0])
    for bc in occ_rates[1:]:
        octxt += ' ' + str('%0.5e' % bc)
    #print(octxt.split()[0]
    return octxt

# makes sure strike angle is within acceptable range
def check_stk_angle(strike):
    
    strike = round(strike)
    if strike >= 360.:
        strike -= 360.
    elif strike < 0.:
        strike += 360
        
    return abs(strike)

# gets nodal plane text for rupture logic tree
def get_nodal_plane_text(m):
    nptxt = ''
    
    # use shmax info assuming shallow dipping reverse faults
    if m['pref_stk'] == -999.0:
        
        # get preferred strikes from shmax for reverse faulting
        pref_stk1 = check_stk_angle(m['shmax']+90)
        pref_stk2 = check_stk_angle(m['shmax']+270)
        
        nptxt += '                <nodalPlane probability="0.34" strike="'+str(pref_stk1)+'" dip="35.0" rake="90.0" />\n'
        nptxt += '                <nodalPlane probability="0.34" strike="'+str(pref_stk2)+'" dip="35.0" rake="90.0" />\n'
        
        # check strike values
        stk_pos1 = check_stk_angle(round(pref_stk1 + m['shmax_sig']))
        stk_neg1 = check_stk_angle(round(pref_stk1 - m['shmax_sig']))
         
        # set shmax +/- sigma
        nptxt += '                <nodalPlane probability="0.08" strike="'+str(stk_pos1)+'" dip="35.0" rake="90.0" />\n'
        nptxt += '                <nodalPlane probability="0.08" strike="'+str(stk_neg1)+'" dip="35.0" rake="90.0" />\n'
        
        # check strike values
        stk_pos2 = check_stk_angle(round(pref_stk2 + m['shmax_sig']))
        stk_neg2 = check_stk_angle(round(pref_stk2 - m['shmax_sig']))
         
        # set shmax +/- sigma
        nptxt += '                <nodalPlane probability="0.08" strike="'+str(stk_pos2)+'" dip="35.0" rake="90.0" />\n'
        nptxt += '                <nodalPlane probability="0.08" strike="'+str(stk_neg2)+'" dip="35.0" rake="90.0" />\n'
            
    # if stk/dip/rke set mostly for offshore regions
    else:
        nptxt += '                <nodalPlane probability="1.0" strike="'+str(m['pref_stk']) \
                                   +'" dip="'+str(m['pref_dip'])+'" rake="'+str(m['pref_rke'])+'" />\n'

    return nptxt

def test_beta_curves(m, binwid, meta, mx_dict):
    from tools.oq_tools import get_oq_incrementalMFD
    import matplotlib.pyplot as plt
    
    mx_vals = mx_dict[m['trt']]['mx_vals']
    
    bbc, bmrange = get_oq_incrementalMFD(m['src_beta'][0], m['src_N0'][0], \
                                        m['min_mag'], mx_vals[2], \
                                        binwid)
    lbc, lmrange = get_oq_incrementalMFD(m['src_beta'][1], m['src_N0'][1], \
                                        m['min_mag'], mx_vals[0], \
                                        binwid)
    ubc, umrange = get_oq_incrementalMFD(m['src_beta'][2], m['src_N0'][2], \
                                        m['min_mag'], mx_vals[4], \
                                        binwid)
    
    plt.semilogy(bmrange, bbc, 'r-')
    plt.semilogy(lmrange, lbc, 'b-')
    plt.semilogy(umrange, ubc, 'g-')
    
    plt.show()

'''
start main code here
'''
def write_oq_sourcefile(model, meta, mx_dict):
    """
    model = a list of dictionaries for each area source
    modelpath = folder for sources to be included in source_model_logic_tree.xml
    logicpath = folder for logic tree
    multimods = argv[2] # for setting weights of alternative models (True or False)
    meta = True gives weight of 1 to best Mmax and b-value
    """

    from oq_tools import beta2bval, get_line_parallels
    from numpy import log10, max, min, tan, radians, isinf, floor
    from os import path
    
    # set big bbox params
    bbmaxlon = -180
    bbmaxlat = -90
    bbminlon = 180
    bbminlat = 90
    
    # make xml header
    header = '<?xml version="1.0" encoding="utf-8"?>\n'
    header += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
    header += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
    
    '''
    # set wieghts
    bval_wt    = [0.68, 0.16, 0.16]
    max_mag_wt = [0.60, 0.30, 0.10]
    '''
    
    outbase = path.split(meta['modelPath'])[-1]
    
    # start xml text
    newxml = header + '    <sourceModel name="'+outbase+'_collapsed">\n\n'
    
    # get src codes and rename if duplicated
    codes = []
    for m in model:
        codes.append(m['src_code'])
    #ucodes = unique(codes)
    
    # start loop thru area sources
    for m in model:
        
        # set magScaleRel
        if float(m['class']) <= 7.:
            magScaleRel = 'Leonard2014_SCR'
            ruptAspectRatio = 1.5 # balance between L14 and Cea14 surface rupture lengths
            min_mag = 4.5
        elif float(m['class']) == 8 or float(m['class']) == 9:
            magScaleRel = 'WC1994'
            ruptAspectRatio = 1.5
            min_mag = 5.5
        elif float(m['class']) == 10:
            magScaleRel = 'StrasserInterface'
            ruptAspectRatio = 1.5 # based on approx AH interface apect ratios at Mw 8
            min_mag = 6.5
        elif floor(float(m['class'])) == 11:
            magScaleRel = 'StrasserIntraslab'
            ruptAspectRatio = 1.2 # based on approx AH intraslab apect ratios at Mw 7.5
            min_mag = 5.5
        
        # comment out sources with null activitiy rates
        if m['src_N0'][-1] == -99.0:
            newxml += '        <!--\n'
        
        #######################################################################
        # write area sources
        #######################################################################
        if m['src_type'] == 'area':
            #print(m['src_type']
            
            # rename source code if "." exists
            m['src_code'].replace('.', '')
            
            newxml += '        <areaSource id="'+m['src_code']+'" name="'+\
                       m['src_name']+'" tectonicRegion="'+m['gmm_trt']+'">\n'
            
            newxml += '            <areaGeometry>\n'
            newxml += '                <gml:Polygon>\n'
            newxml += '                    <gml:exterior>\n'
            newxml += '                        <gml:LinearRing>\n'
            newxml += '                            <gml:posList>\n'
                            
            # get polygon text
            polytxt = ''
            pp = 0
            for xy in m['src_shape'][:-1]: # no need to close poly
                addPoint = True
                # check if duplicating points
                if pp > 0:
                   if xy[0] == xy0[0] and xy[1] == xy0[1]:
                      addPoint = False
                      
                if addPoint == True:
                    polytxt = polytxt + '                                ' + str("%0.4f" % xy[0]) \
                                      + ' ' + str("%0.4f" % xy[1]) + '\n'
                
                xy0 = xy
                pp += 1
                
            # add poly text
            newxml += polytxt
            
            newxml += '                            </gml:posList>\n'
            newxml += '                        </gml:LinearRing>\n'
            newxml += '                    </gml:exterior>\n'
            newxml += '                </gml:Polygon>\n'
            
            ###################################################################
            # print(model bbox of model
            
            # this is not required for the nrml files, but useful for setting up job.ini files
            
            buff = 0.1
            maxlon = max(m['src_shape'][:,0])+buff
            minlon = min(m['src_shape'][:,0])-buff
            maxlat = max(m['src_shape'][:,1])+buff
            minlat = min(m['src_shape'][:,1])-buff
    
            # get big bbox
            if maxlon > bbmaxlon: bbmaxlon = maxlon
            if minlon < bbminlon: bbminlon = minlon
            if maxlat > bbmaxlat: bbmaxlat = maxlat
            if minlat < bbminlat: bbminlat = minlat
    
            #print(m['src_code'], minlon, minlat, ',', minlon, maxlat, ',', maxlon, maxlat, ',', maxlon, minlat
            ###################################################################
    
            # set depth distribution
            if m['src_dep'][0] <= m['src_usd'] or m['src_dep'][0] >= m['src_lsd']:
                print(m['src_code'], 'FIX DEPTHS')
                
            newxml += '                <upperSeismoDepth>'+str(m['src_usd'])+'</upperSeismoDepth>\n'
            newxml += '                <lowerSeismoDepth>'+str(m['src_lsd'])+'</lowerSeismoDepth>\n'
            
            # set source geometry
            newxml += '            </areaGeometry>\n'
            newxml += '            <magScaleRel>'+magScaleRel+'</magScaleRel>\n'
            newxml += '            <ruptAspectRatio>'+str(ruptAspectRatio)+'</ruptAspectRatio>\n'
            
            # get weighted rates
            binwid = 0.1
            octxt = make_collapse_occurrence_text(m, min_mag, binwid, meta, mx_dict)
                                 
            newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (min_mag+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
            newxml += '                <occurRates>'+octxt+'</occurRates>\n'
            newxml += '            </incrementalMFD>\n'
            
            """
            # set GR recurrence pars
            tmpN0   = m['src_N0'][i]
            tmpbeta = m['src_beta'][i]
            tmpmmax = m['max_mag'][j]
            grtxt = ''.join(('            <truncGutenbergRichterMFD aValue="', \
                            str("%0.4f" % log10(tmpN0)),'" bValue="', \
                            str("%0.4f" % beta2bval(tmpbeta)),'" minMag="', \
                            str("%0.2f" % m['min_mag']),'" maxMag="', \
                            str("%0.2f" % tmpmmax),'"/>\n'))
                                
            newxml += grtxt
            """
            # set nodal planes
            newxml += '            <nodalPlaneDist>\n'
            newxml += get_nodal_plane_text(m)
            newxml += '            </nodalPlaneDist>\n'
            
            # set hypo depth
            newxml += '            <hypoDepthDist>\n'
            if m['src_dep'][1] != -999.0:
                newxml += '                <hypoDepth probability="0.50" depth="'+str("%0.1f" % m['src_dep'][0])+'"/>\n' \
                         +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][1])+'"/>\n' \
                         +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][2])+'"/>\n'
            else:
                newxml += '                <hypoDepth probability="1.0" depth="'+str("%0.1f" % m['src_dep'][0])+'"/>\n'
            
            newxml += '            </hypoDepthDist>\n'
            
            if m['src_N0'][-1] == -99.0:
                newxml += '        </areaSource>\n'
            else:
                newxml += '        </areaSource>\n\n'
        #######################################################################
        # now make fault sources
        #######################################################################
        elif m['src_type'] == 'fault':
             
            # rename source code if "." exists
            m['src_code'].replace('.', '')
            src_code = m['src_code']
            
            if isinf(log10(m['src_N0'][0])) == False:
                ###################################################################
                # do complex faults
                ###################################################################
                if m['fault_dip'][0] != m['fault_dip'][1]:
                #if m['fault_dip'][0] >= 0: # catches all faults
                    #if m['fault_dip'][0] > 0:
    
                    # id subcript
                    idsub = str("%0.1f" % beta2bval(m['src_beta'][0]))
                    idsub = idsub.replace(".", "")
                    
                    newxml += '        <complexFaultSource id="'+src_code+idsub+'" name="'+\
                               m['src_name']+'" tectonicRegion="'+m['gmm_trt']+'">\n'
                    newxml += '            <complexFaultGeometry>\n'
                    newxml += '                <faultTopEdge>\n'
                    newxml += '                    <gml:LineString>\n'
                    newxml += '                        <gml:posList>\n'
                        
                    # calculate lat lons from surface projection
                    # get upper h-dist
                    upperhdist = m['src_dep'][0] / tan(radians(m['fault_dip'][0]))
                    upperxy = get_line_parallels(m['src_shape'], upperhdist)[0]
                
                    # make upper text
                    xytxt = ''
                    for xy in upperxy:
                        xytxt += '                            ' + \
                                 ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][0])))+'\n'
                    newxml += xytxt
                    newxml += '                        </gml:posList>\n'
                    newxml += '                    </gml:LineString>\n'
                    newxml += '                </faultTopEdge>\n'
                    newxml += '                <intermediateEdge>\n'
                    newxml += '                    <gml:LineString>\n'
                    newxml += '                        <gml:posList>\n'
                
                    
                    # calculate lat lons from upper edge
                    # get intermediate h-dist
                    interhdist = (m['src_dep'][1] - m['src_dep'][0]) / tan(radians(m['fault_dip'][0]))
                    interxy = get_line_parallels(upperxy, interhdist)[0]
                
                    # make intermediate text
                    xytxt = ''
                    for xy in interxy:
                        xytxt += '                            ' + \
                                 ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][1])))+'\n'
                    newxml += xytxt
                    newxml += '                        </gml:posList>\n'
                    newxml += '                    </gml:LineString>\n'
                    newxml += '                </intermediateEdge>\n'
                    newxml += '                <faultBottomEdge>\n'
                    newxml += '                    <gml:LineString>\n'
                    newxml += '                        <gml:posList>\n'
                
                    # calculate lat lons from intermediate edge
                    # get bottom h-dist
                    bottomhdist = (m['src_dep'][2] - m['src_dep'][1]) / tan(radians(m['fault_dip'][1]))
                    bottomxy = get_line_parallels(interxy, bottomhdist)[0]
                
                    # make bottom text
                    xytxt = ''
                    for xy in bottomxy:
                        xytxt += '                            ' + \
                                 ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1]), str(m['src_dep'][2])))+'\n'
                    newxml += xytxt
                    newxml += '                        </gml:posList>\n'
                    newxml += '                    </gml:LineString>\n'
                    newxml += '                </faultBottomEdge>\n'
                    newxml += '            </complexFaultGeometry>\n'
                    
                    '''
                    # get fault area scaling model
                    '''
                    #src_code = m['src_code']
                    if src_code.startswith('CIS'):
                        newxml += '            <magScaleRel>GSCCascadia</magScaleRel>\n'
                    elif src_code.startswith('WIN'):
                        newxml += '            <magScaleRel>GSCOffshoreThrustsWIN</magScaleRel>\n'
                    elif src_code.startswith('HGT'):
                        newxml += '            <magScaleRel>GSCOffshoreThrustsHGT</magScaleRel>\n'
                    elif src_code.startswith('QCSS') or src_code.startswith('FWF'):
                        newxml += '            <magScaleRel>WC1994_QCSS</magScaleRel>\n'
                    elif src_code.startswith('EISO'):
                        newxml += '            <magScaleRel>GSCEISO</magScaleRel>\n'
                    elif src_code.startswith('EISB'):
                        newxml += '            <magScaleRel>GSCEISB</magScaleRel>\n'
                    elif src_code.startswith('EISI'):
                        newxml += '            <magScaleRel>GSCEISI</magScaleRel>\n'
                    else:
                        newxml += '            <magScaleRel>'+magScaleRel+'</magScaleRel>\n'
                        #newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
                    
                    newxml += '            <ruptAspectRatio>'+str(ruptAspectRatio)+'</ruptAspectRatio>\n'
                
                    '''
                    # now get appropriate MFD
                    '''
                    # do incremental MFD
                    if m['src_beta'][0] > -99:
                        # adjust N0 value to account for weighting of fault sources
                    
                        octxt = make_collapse_occurrence_text(m, min_mag, binwid, meta, mx_dict)
                                    
                        # make text
                        newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                        newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                        newxml += '            </incrementalMFD>\n'
                
                    
                    if m['fault_dip'][0] != 90.:
                        newxml += '            <rake>90.0</rake>\n'
                    else:
                        newxml += '            <rake>0.0</rake>\n'
                
                    newxml += '        </complexFaultSource>\n\n'
    
                ###################################################################
                # else do simple fault
                ###################################################################
                elif m['fault_dip'][0] == m['fault_dip'][1]:
                    
                    # id subcript
                    idsub = str("%0.1f" % beta2bval(m['src_beta'][0]))
                    idsub = idsub.replace(".", "")
                    
                    newxml += '        <simpleFaultSource id="'+m['src_code']+idsub+'" name="'+\
                                         m['src_name']+'" tectonicRegion="'+m['gmm_trt']+'">\n'
                    newxml += '            <simpleFaultGeometry>\n'
                    newxml += '                <gml:LineString>\n'
                    newxml += '                    <gml:posList>\n'
                
                    # simple fauls use surface projection!
                    '''
                    # calculate lat lons from surface projection
                    # get upper h-dist
                    upperhdist = m['src_dep'][0] / tan(radians(m['fault_dip'][0]))
                    upperxy = get_line_parallels(m['src_shape'], upperhdist)[0]
                    '''
                    
                    xytxt = ''
                    for xy in m['src_shape']:
                        xytxt += '                            ' + \
                                 ' '.join((str('%0.4f' % xy[0]), str('%0.4f' % xy[1])))+'\n'
                    newxml += xytxt
                    
                    newxml += '                    </gml:posList>\n'
                    newxml += '                </gml:LineString>\n'
                    newxml += '                <dip>'+str(m['fault_dip'][0])+'</dip>\n'
                    newxml += '                <upperSeismoDepth>'+str(m['src_dep'][0])+'</upperSeismoDepth>\n'
                    newxml += '                <lowerSeismoDepth>'+str(m['src_dep'][-1])+'</lowerSeismoDepth>\n'
                    newxml += '            </simpleFaultGeometry>\n'
                    
                    '''
                    # get fault area scaling model
                    '''
                    src_code = m['src_code']
                    if src_code == 'CIS':
                        newxml += '            <magScaleRel>GSCCascadia</magScaleRel>\n'
                    elif src_code.startswith('WIN'):
                        newxml += '            <magScaleRel>GSCOffshoreThrustsWIN</magScaleRel>\n'
                    elif src_code.startswith('HGT'):
                        newxml += '            <magScaleRel>GSCOffshoreThrustsHGT</magScaleRel>\n'
                    elif src_code.startswith('QCSS') or src_code.startswith('FWF'):
                        newxml += '            <magScaleRel>WC1994_QCSS</magScaleRel>\n'
                    elif src_code.startswith('EISO'):
                        newxml += '            <magScaleRel>GSCEISO</magScaleRel>\n'
                    elif src_code.startswith('EISB'):
                        newxml += '            <magScaleRel>GSCEISB</magScaleRel>\n'
                    elif src_code.startswith('EISI'):
                        newxml += '            <magScaleRel>GSCEISI</magScaleRel>\n'
                    else:
                        newxml += '            <magScaleRel>'+magScaleRel+'</magScaleRel>\n'
                        
                    newxml += '            <ruptAspectRatio>'+str(ruptAspectRatio)+'</ruptAspectRatio>\n'
                    
                    '''
                    # now get appropriate MFD
                    '''
                    # do incremental MFD
                    if m['src_beta'][0] > -99:
                        
                        octxt = make_collapse_occurrence_text(m, min_mag, binwid, meta, mx_dict)
                                    
                        # make text
                        newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
                        newxml += '                <occurRates>'+octxt+'</occurRates>\n'
                        newxml += '            </incrementalMFD>\n'
                                    
                    if m['fault_dip'][0] != 90.:
                        newxml += '            <rake>90.0</rake>\n'
                    else:
                        newxml += '            <rake>0.0</rake>\n'
                    
                    if m['src_N0'][-1] == -99.0:
                        newxml += '        </simpleFaultSource>\n'
                    else:
                        newxml += '        </simpleFaultSource>\n\n'
                
        # comment sources with null activity rates
        if m['src_N0'][-1] == -99.0:
            newxml += '        -->\n\n'
            
    ######################################################################
    # add Australian fault-source model
    ######################################################################
    if meta['doSeisTec'] == True:
        aust_fault_file = path.join('..', 'faults', 'National_Fault_Source_Model_2018_Collapsed_NSHA13', \
                                    'National_Fault_Source_Model_2018_Collapsed_NSHA13_all_methods_collapsed_inc_cluster_gmm_trt.xml')
        lines = open(aust_fault_file).readlines()[3:-2]
        for line in lines:
            newxml += '    ' + line
    
    ######################################################################
    # add indoneasia-png area and fault-source model
    ######################################################################
    
    #indo_png_fault_file = path.join('..', 'banda', 'Banda_Fault_Sources_NSHA_2018.xml')
    indo_png_source_file = path.join('2018_mw', 'Java_Banda_PNG', 'input', 'collapsed', 'Java_Banda_PNG_collapsed.xml')
    lines = open(indo_png_source_file).readlines()[4:-2]
    for line in lines:
        newxml += line
    
    #print('\nSkipping Banda Faults\n'
       
    ######################################################################
    # finish nrml
    newxml += '    </sourceModel>\n'
    newxml += '</nrml>'
    
    # write Big BBOX
    #print('\nBBOX:', bbminlon, bbminlat, ',', bbminlon, bbmaxlat, ',', bbmaxlon, bbmaxlat, ',', bbmaxlon, bbminlat
    
    # write new data to file
    outxml = path.join(meta['modelPath'], meta['modelFile'])
    #outxml = '/'.join((src_folder, ''.join((outbase,'_',bl,'_',ml,'.xml'))))
    f = open(outxml,'w')
    f.write(newxml)
    f.close()
    
    #srcxmls.append(outxml)
    return outxml

######################################################################
# now that the source file have been written, make the logic tree file
######################################################################

def make_logic_tree(srcxmls, branch_wts, meta):    
    from os import path
    
    #print(branch_wts
    # if multimodel - adjust weights
    '''
    if meta['multiMods'] == 'True':
        branch_wts = array(branch_wts)
        branch_wt *= m['src_reg_wt']
        print('Branch Weights: ', m['src_reg_wt']
        #else:
        #    full_wt = concatenate((branch_wt, branch_wt, branch_wt))
    '''
    
    newxml = '<?xml version="1.0" encoding="UTF-8"?>\n'
    newxml += '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
    newxml += '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n\n'
    newxml += '    <logicTree logicTreeID="lt1">\n'
    newxml += '        <logicTreeBranchingLevel branchingLevelID="bl1">\n'
    newxml += '            <logicTreeBranchSet uncertaintyType="sourceModel"\n' \
          '                                branchSetID="bs1">\n\n'
    
    # make branches
    for i, srcxml in enumerate(srcxmls):
        #print(i, srcxml
        #logictreepath = logicpath + sep + path.split(branch)[-1]
        if meta['splitXMLPath'] == True:
            logictreepath = path.split(srcxml)[-1]
        else:
            logictreepath = srcxml
        #print(i, logictreepath
            
        newxml += '                <logicTreeBranch branchID="b' + str(i+1) + '">\n'
        newxml += '                    <uncertaintyModel>'+logictreepath+'</uncertaintyModel>\n'
        newxml += '                    <uncertaintyWeight>'+str(branch_wts[i])+'</uncertaintyWeight>\n'
        newxml += '                </logicTreeBranch>\n\n'
    
    newxml += '            </logicTreeBranchSet>\n'
    newxml += '        </logicTreeBranchingLevel>\n'
    newxml += '    </logicTree>\n'
    newxml += '</nrml>'
        
    # write logic tree to file
    outxml = path.join(meta['modelPath'], 
                       ''.join((meta['modelFile'].split('.')[0].split('_')[0], #[:-10], 
                                '_source_model_logic_tree.xml')))
                                
    f = open(outxml,'w')
    f.write(newxml)
    f.close()
    
    
###############################################################################
# make source dict for OQ input writer
###############################################################################
# set code prefix to ensure unique ids and optimise jobs
def get_code_prefix(modelshp):
    print(modelshp)
    if modelshp.endswith('ARUP_NSHA18_MFD.shp') or modelshp.endswith('ARUP_NSHA18_MFD_MX.shp'):
        code_prefix = 'ARUP'
    elif modelshp.endswith('ARUP_Background_NSHA18_MFD.shp') or modelshp.endswith('ARUP_Background_NSHA18_MFD_MX.shp'):
        code_prefix = 'ARUPB'
    elif modelshp.endswith('AUS6_NSHA18_MFD.shp') or modelshp.endswith('AUS6_NSHA18_MFD_MX.shp'):
        code_prefix = 'AUS6'
    elif modelshp.endswith('DIMAUS_NSHA18_MFD.shp') or modelshp.endswith('DIMAUS_NSHA18_MFD_MX.shp'):
        code_prefix = 'DIM'
    elif modelshp.endswith('Domains_NSHA18_MFD.shp') or modelshp.endswith('Domains_NSHA18_MFD_MX.shp'):
        code_prefix = 'DOM'
    elif modelshp.endswith('Leonard2008_NSHA18_MFD.shp') or modelshp.endswith('Leonard2008_NSHA18_MFD_MX.shp'):
        code_prefix = 'L08'
    elif modelshp.endswith('NSHA13_NSHA18_MFD.shp') or modelshp.endswith('NSHA13_NSHA18_MFD_MX.shp'):
        code_prefix = 'NSHM'
    elif modelshp.endswith('NSHA13_BACKGROUND_NSHA18_MFD.shp') or modelshp.endswith('NSHA13_BACKGROUND_NSHA18_MFD_MX.shp'):
        code_prefix = 'NSHMB'
    elif modelshp.endswith('SIN_MCC_NSHA18_MFD.shp') or modelshp.endswith('SIN_MCC_NSHA18_MFD_MX.shp'):
        code_prefix = 'SM'
    elif modelshp.endswith('AUS6_Gridded_b_MFD.shp'):
        code_prefix = 'AUS6G'
    elif modelshp.endswith('DIMAUS_Gridded_b_MFD.shp'):
        code_prefix = 'DIMG'
    elif modelshp.endswith('NSHA13_Gridded_b_MFD.shp'):
        code_prefix = 'NSHMG'
    else:
        code_prefix = ''
        
    return code_prefix
    
def src_shape2dict(modelshp):
    import shapefile  
    from numpy import array
    #from oq_tools import get_oq_incrementalMFD
    from tools.nsha_tools import bval2beta
    
    # read new shapefile
    sf = shapefile.Reader(modelshp)
    records = sf.records()
    shapes  = sf.shapes()

    # set model list
    model = []
    
    # set source code prefix
    code_prefix = get_code_prefix(modelshp)
    
    # loop thru recs and make dict
    for rec, shape in zip(records, shapes):
        if not float(rec[15]) == -99:
            # determine whether source optimization occurs - ignores banda, etc sources to ensure same params are used
            if float(rec[3]) <= 8.:
                optim_code = '_'.join((code_prefix, rec[1]))
                
                m = {'src_name':rec[0], 'src_code':optim_code, 'src_type':rec[2],
                     'class':rec[3], 'trt':rec[33], 'src_shape':array(shape.points),
                     'src_dep':[float(rec[6]), float(rec[7]), float(rec[8])],
                     'src_usd':float(rec[9]), 'src_lsd':float(rec[10]),
                     'max_mag':[float(rec[14]), float(rec[15]), float(rec[16])],
                     'src_N0':[float(rec[17]), float(rec[18]), float(rec[19])],
                     'src_beta':[bval2beta(float(rec[20])), bval2beta(float(rec[21])), bval2beta(float(rec[22]))],
                     'min_mag':float(rec[12]), 'src_weight':float(rec[4]), 'src_reg_wt':1.,
                     'rate_adj_fact':float(rec[5]), 'pref_stk':float(rec[28]), 'pref_dip':float(rec[29]),
                     'pref_rke':float(rec[30]), 'shmax':float(rec[31]), 'shmax_sig':float(rec[32]), 'gmm_trt':rec[34]}
                
                model.append(m)
                
            '''
            else:
                optim_code = rec[1]
            '''    
            
            
            
        else:
            # parse GSC version
            m = {'src_name':rec[0], 'src_code':rec[1], 'src_type':rec[2], 'src_type':rec[2],
                  'src_weight':float(rec[3]), 'src_shape':array(shape.points), 
                  'src_dep':[float(rec[4]), float(rec[5]), float(rec[6])], 'min_mag':float(rec[7]), 
                  'max_mag':[float(rec[9]), float(rec[10]), float(rec[11])],
                  'src_N0':[float(rec[12]), float(rec[13]), float(rec[14])],
                  'src_beta':[float(rec[15]), float(rec[16]), float(rec[17])],
                  'src_reg_wt':1, 'trt':rec[23]}
    
    return model 
"""    
0     w.field('SRC_NAME','C','100')
1     w.field('CODE','C','10')
2     w.field('SRC_TYPE','C','10')
3     w.field('CLASS','C','10')
4     w.field('SRC_WEIGHT','F', 8, 2)
5     w.field('RTE_ADJ_F','F', 6, 4)
6     w.field('DEP_BEST','F', 6, 1)
7     w.field('DEP_UPPER','F', 6, 1)
8     w.field('DEP_LOWER','F', 6, 1)
9     w.field('USD','F', 4, 1)
10    w.field('LSD','F', 4, 1)
11    w.field('OW_LSD','F', 4, 1)
12    w.field('MIN_MAG','F', 4, 2)
13    w.field('MIN_RMAG','F', 4, 2)
14    w.field('MMAX_BEST','F', 4, 2)
15    w.field('MMAX_LOWER','F', 4, 2)
16    w.field('MMAX_UPPER','F', 4, 2)
17    w.field('N0_BEST','F', 8, 5)
18    w.field('N0_LOWER','F', 8, 5)
19    w.field('N0_UPPER','F', 8, 5)
20    w.field('BVAL_BEST','F', 6, 3)
21    w.field('BVAL_LOWER','F', 6, 3)
22    w.field('BVAL_UPPER','F', 6, 3)
23    w.field('BVAL_FIX','F', 6, 3)
24    w.field('BVAL_FIX_S','F', 6, 3)
25    w.field('YCOMP','C','70')
26    w.field('MCOMP','C','50')
27    w.field('CAT_YMAX', 'F', 8, 3)
28    w.field('PREF_STK','F', 6, 2)
29    w.field('PREF_DIP','F', 6, 2)
30    w.field('PREF_RKE','F', 6, 2)
31    w.field('SHMAX','F', 6, 2)
32    w.field('SHMAX_SIG','F', 6, 2)
33    w.field('TRT','C','100')
34    w.field('GMM_TRT','C','100')
35    w.field('DOMAIN','F', 2, 0)
36    w.field('CAT_FILE','C','50')
"""