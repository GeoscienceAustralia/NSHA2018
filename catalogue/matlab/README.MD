Steps to generate matlab mat file, to be used in magnitude revision and conversion routines

1) Parse Gary Gibson's earthquake catalogue using "ggcat2mat.m"
2) Merge GA magnitudes for larger earthquakes using "mergeGA.m"
	- semi-automated code, requires user input when multiple events found
3) Append alternate ML values using "merge_alt_mags.m"
4) Append update MW values using "append_mw.m"
	- exports "mdat_mw.mat"
5) Get preferred ML and MW magnitudes by running "get_pref_mags.m"
	- assumes equivalence between ML and MD, MP and M?
	- exports "mdat_pref_mag_types.mat"
6) Make consistent magnitudes across all time periods using "reviseMLs_2018.m"
	- exports "mdat_ml_rev.mat"
7) Make generic ML correction factors for those events that don't get caught by the logic in Step 6 by running "get_M_corr_fact.m"
	- This should only need to be run once
8) Re-run "reviseMLs_2018.m" to propagate any changes from Step 7 to catalogue
9) Make preferred MW catalogue using "get_pref_MW.m"
	- uses polynomial ML2MW conversion from simulated data
	- exports "mdat_mw_pref.mat"
