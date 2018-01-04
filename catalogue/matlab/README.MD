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
6) Make consistent magnitudes across all time periods using 'reviseMLs_2018.m"
	- exports "mdat_ml_rev.mat"
