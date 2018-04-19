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