#####################################################################################
# imports source model from get_NSHA18_mfds.py and outputs an OpenQuake source file with
# collapsed earthquake rates

# Usage: python make_openquake_source_file.py <path to pkl file> <output base name>
#            e.g. python make_openquake_source_file.py SWCan_T3EclC1.pkl swcan
#####################################################################################

def make_collapse_occurrence_text(m, binwid, meta, mx_dict):
    from numpy import array, zeros, argmax, zeros_like
    from oq_tools import get_oq_incrementalMFD
    '''
    best_beta: boolean; True if only using best curve
    mx_vals: array of Mmax values
    mx_wts: array of Mmax weights
    '''
    
    # set mx arrays from mx dictionary using zone trt
    mx_vals = mx_dict[m['trt']]['mx_vals']
    mx_wts  = mx_dict[m['trt']]['mx_wts']
    
    # REMOVE ME WHEN FINISHED TESTING
    mx_wts = [0.1, 0.2, 0.4, 0.2, 0.1]
        
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
        effN0 = N0 * m['src_weight']
        
        for mx, mxwt in zip(mx_vals, mx_wts):
            
            betacurve, mrange = get_oq_incrementalMFD(beta_val, effN0, \
                                                      m['min_mag'], mx, \
                                                      binwid)
            
            #print beta_val, beta_wt, N0, mx, mxwt, beta_wt*mxwt
            
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
                
    # convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(wtd_rates[0:-1])):
        occ_rates.append(wtd_rates[b] - wtd_rates[b+1])
    occ_rates.append(wtd_rates[-1])
    
    # make text object                        
    octxt = str('%0.5e' % occ_rates[0])
    for bc in occ_rates[1:]:
        octxt += ' ' + str('%0.5e' % bc)
    
    print occ_rates[0]
    return octxt

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
    from numpy import log10, max, min, tan, radians, isinf
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
        # comment out sources with null activitiy rates
        if m['src_N0'][-1] == -99.0:
            newxml += '        <!--\n'
        
        #######################################################################
        # write area sources
        #######################################################################
        if m['src_type'] == 'area':
            print m['src_type']
            
            # rename source code if "." exists
            m['src_code'].replace('.', '')
            
            newxml += '        <areaSource id="'+m['src_code']+'" name="'+\
                       m['src_name']+'" tectonicRegion="'+m['trt']+'">\n'
            
            newxml += '            <areaGeometry>\n'
            newxml += '                <gml:Polygon>\n'
            newxml += '                    <gml:exterior>\n'
            newxml += '                        <gml:LinearRing>\n'
            newxml += '                            <gml:posList>\n'
                            
            # get polygon text
            polytxt = ''
            for xy in m['src_shape'][:-1]: # no need to close poly
                polytxt = polytxt + '                                ' + str("%0.4f" % xy[0]) \
                                  + ' ' + str("%0.4f" % xy[1]) + '\n'
            newxml += polytxt
            
            newxml += '                            </gml:posList>\n'
            newxml += '                        </gml:LinearRing>\n'
            newxml += '                    </gml:exterior>\n'
            newxml += '                </gml:Polygon>\n'
            
            ###################################################################
            # print model bbox of model
            
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
    
            print m['src_code'], minlon, minlat, ',', minlon, maxlat, ',', maxlon, maxlat, ',', maxlon, minlat
            ###################################################################
    
            # set depth distribution
            if min(m['src_dep']) != max(m['src_dep']):
                #newxml += '                <upperSeismoDepth>'+str("%0.1f" % min(m['src_dep']))+'</upperSeismoDepth>\n'
                #newxml += '                <lowerSeismoDepth>'+str("%0.1f" % max(m['src_dep']))+'</lowerSeismoDepth>\n'
                newxml += '                <upperSeismoDepth>0.0</upperSeismoDepth>\n'
                newxml += '                <lowerSeismoDepth>20.0</lowerSeismoDepth>\n'
            else:
                newxml += '                <upperSeismoDepth>'+str("%0.1f" % (min(m['src_dep'])-10))+'</upperSeismoDepth>\n'
                newxml += '                <lowerSeismoDepth>'+str("%0.1f" % (min(m['src_dep'])+10))+'</lowerSeismoDepth>\n'
                
            newxml += '            </areaGeometry>\n'
            newxml += '            <magScaleRel>Leonard2014_SCR</magScaleRel>\n'
            #newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
            #newxml += '            <ruptAspectRatio>2.0</ruptAspectRatio>\n'
            newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
            
            # get weighted rates
            binwid = 0.1
            octxt = make_collapse_occurrence_text(m, binwid, meta, mx_dict)
                                 
            newxml += '            <incrementalMFD minMag="'+str('%0.2f' % (m['min_mag']+0.5*binwid))+'" binWidth="'+str(binwid)+'">\n'
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
            
            newxml += '                <nodalPlane probability="0.3" strike="0.0" dip="30.0" rake="90.0" />\n'
#            newxml += '                <nodalPlane probability="0.0625" strike="45.0" dip="30.0" rake="90.0" />\n'
            newxml += '                <nodalPlane probability="0.2" strike="90.0" dip="30.0" rake="90.0" />\n'
#            newxml += '                <nodalPlane probability="0.0625" strike="135.0" dip="30.0" rake="90.0" />\n'
            newxml += '                <nodalPlane probability="0.3" strike="180.0" dip="30.0" rake="90.0" />\n'
#            newxml += '                <nodalPlane probability="0.0625" strike="225.0" dip="30.0" rake="90.0" />\n'
            newxml += '                <nodalPlane probability="0.2" strike="270.0" dip="30.0" rake="90.0" />\n'
#            newxml += '                <nodalPlane probability="0.0625" strike="315.0" dip="30.0" rake="90.0" />\n'
    
            newxml += '            </nodalPlaneDist>\n'
            
    
            # set hypo depth
            newxml += '            <hypoDepthDist>\n'
            newxml += '                <hypoDepth probability="0.50" depth="'+str("%0.1f" % m['src_dep'][0])+'"/>\n' \
                     +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][1])+'"/>\n' \
                     +'                <hypoDepth probability="0.25" depth="'+str("%0.1f" % m['src_dep'][2])+'"/>\n'
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
                               m['src_name']+'" tectonicRegion="'+m['trt']+'">\n'
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
                        newxml += '            <magScaleRel>Leonard2014_SCR</magScaleRel>\n'
                        #newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
                    
                    newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
                
                    '''
                    # now get appropriate MFD
                    '''
                    # do incremental MFD
                    if m['src_beta'][0] > -99:
                        # adjust N0 value to account for weighting of fault sources
                    
                        octxt = make_collapse_occurrence_text(m, binwid, meta, mx_dict)
                                    
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
                                         m['src_name']+'" tectonicRegion="'+m['trt']+'">\n'
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
                        newxml += '            <magScaleRel>Leonard2014_SCR</magScaleRel>\n'
                        #newxml += '            <magScaleRel>WC1994</magScaleRel>\n'
                    
                    newxml += '            <ruptAspectRatio>1.0</ruptAspectRatio>\n'
                    #newxml += '            <ruptAspectRatio>2.0</ruptAspectRatio>\n'
                    '''
                    # now get appropriate MFD
                    '''
                    # do incremental MFD
                    if m['src_beta'][0] > -99:
                        
                        octxt = make_collapse_occurrence_text(m, binwid, meta, mx_dict)
                                    
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
    
    # finish nrml
    newxml += '    </sourceModel>\n'
    newxml += '</nrml>'
    
    # write Big BBOX
    print '\nBBOX:', bbminlon, bbminlat, ',', bbminlon, bbmaxlat, ',', bbmaxlon, bbmaxlat, ',', bbmaxlon, bbminlat
    
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
    
    # if multimodel - adjust weights
    '''
    if meta['multiMods'] == 'True':
        branch_wts = array(branch_wts)
        branch_wt *= m['src_reg_wt']
        print 'Branch Weights: ', m['src_reg_wt']
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
        #logictreepath = logicpath + sep + path.split(branch)[-1]
        logictreepath = path.split(srcxml)[-1]
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
                       ''.join((meta['modelFile'].split('.')[0][:-11], 
                                '_source_model_logic_tree.xml')))
                                
    f = open(outxml,'w')
    f.write(newxml)
    f.close()
    
    
###############################################################################
# make source dict for OQ input writer
###############################################################################

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
    
    # loop thru recs and make dict
    for rec, shape in zip(records, shapes):
        if not float(rec[15]) == -99:
            m = {'src_name':rec[0], 'src_code':rec[1], 'src_type':rec[2], 
                 'trt':rec[23], 'src_shape':array(shape.points), 
                 'src_dep':[float(rec[4]), float(rec[5]), float(rec[6])], 
                 'src_N0':[float(rec[12]), float(rec[13]), float(rec[14])], 
                 'src_beta':[bval2beta(float(rec[15])), bval2beta(float(rec[16])), bval2beta(float(rec[17]))], 
                 'max_mag':[float(rec[9]), float(rec[10]), float(rec[11])], 
                 'min_mag':float(rec[7]), 'src_weight':float(rec[3]), 'src_reg_wt':1}
            model.append(m)
    
    return model
