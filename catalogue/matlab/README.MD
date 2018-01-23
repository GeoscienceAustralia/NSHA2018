Steps to generate matlab mat file, to be used in magnitude revision and conversion routines

1)  Parse Gary Gibson's earthquake catalogue using "ggcat2mat.m"
2)  Merge GA magnitudes for larger (ML >= 3.25) earthquakes using "mergeGA.m"
	- semi-automated code, requires user input when multiple events found
3)  Merge ISC magnitudes for larger (ML >= 3.25) earthquakes using "mergeISC.m"
	- semi-automated code, requires user input when multiple events found
	- appends USGS PDE post 2015-01-31
4)  Append alternate ML values using "merge_alt_mags.m"
5)  Append MW values using "append_mw.m"
	- exports "mdat_mw.mat"
6)  Get preferred ML and MW magnitudes by running "get_pref_mags.m"
	- assumes equivalence between ML and MD, MP and M?
	- exports "mdat_pref_mag_types.mat"
7)  Make consistent magnitudes across all time periods using "reviseMLs_2018.m"
	- exports "mdat_ml_rev.mat"
8)  Make generic ML correction factors for those events that do not get caught by the logic in Step 6 by running "get_M_corr_fact.m"
	- This should only need to be run once
9)  Re-run "reviseMLs_2018.m" to propagate any changes from Step 7 to catalogue
10) Make preferred MW catalogue using "get_pref_MW.m"
	- uses polynomial ML2MW conversion from simulated data
	- exports "mdat_mw_pref.mat"
	- exports raw catalogue to "../data/NSHA18CAT.MW.VX.X.csv"
11) Make catalogue files for magnitude conversions using "export_mag_comparisons.m":
	- exports NSHA18CAT.mb-MW.csv
	- exports NSHA18CAT.MS-MW.csv
	- exports NSHA18CAT.ML-MW.csv
