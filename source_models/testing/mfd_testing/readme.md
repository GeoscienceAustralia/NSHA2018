For all tests, assume:
	- N0   = 4e+04
	- a    = 4.6020599913279625
	- b    = 1.2
	- beta = 2.7631021115928549

truncGutenbergRichterMFD, binWidth = 0.1 & 0.01
	
<truncGutenbergRichterMFD aValue="4.6020599913279625" bValue="1.2" minMag="4.5" maxMag="7.5" />

For incremental, use

from tools.oq_tools import get_oq_incrementalMFD
betacurve, mrange = get_oq_incrementalMFD(beta_val, effN0, min_mag, mx, binwid)

# convert cummulative rates to annual occurrence rates
    occ_rates = []
    for b in range(0, len(wtd_rates[0:-1])):
        occ_rates.append(wtd_rates[b] - wtd_rates[b+1])
    occ_rates.append(wtd_rates[-1])
    
    # make text object                        
    octxt = str('%0.5e' % occ_rates[0])
    for bc in occ_rates[1:]:
        octxt += ' ' + str('%0.5e' % bc)
    #print octxt.split()[0]
    return octxt

betacurve, mrange = get_oq_incrementalMFD(bval2beta(1.2), 4e+04, 4.5, 7.5, 0.1)

array([1.59202868e-01, 1.20758069e-01, 9.15947062e-02, 6.94720332e-02,
       5.26902696e-02, 3.99600000e-02, 3.03031030e-02, 2.29775975e-02,
       1.74206333e-02, 1.32052449e-02, 1.00075457e-02, 7.58184288e-03,
       5.74175909e-03, 4.34591279e-03, 3.28705509e-03, 2.48382938e-03,
       1.87452037e-03, 1.41231222e-03, 1.06169148e-03, 7.95718453e-04,
       5.93957278e-04, 4.40905774e-04, 3.24804336e-04, 2.36732389e-04,
       1.69922984e-04, 1.19242868e-04, 8.07980689e-05, 5.16347062e-05,
       2.95120332e-05, 1.27302696e-05])


betacurve, mrange = get_oq_incrementalMFD(bval2beta(1.2), 4e+04, 4.5, 7.5, 0.01)

oq-engine --rh job_locs_PGA.ini --exports csv

