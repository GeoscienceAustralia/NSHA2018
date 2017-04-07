The 2012 NSHA catalogue is used as the starting poitnt for the development of the draft design hazard values as prepared for the 2019 AS1170.4.  Geoscience Australia intend to update the catalogue prior to the finalisation of seismic design values in late 2017.  The steps used to develop the catalogue are:

1) Use 2012 Catalogue V0.11 as the starting point (see Burbidge, 2012 Appendix A)
2) Update ML correction code to:
	- replace Richter (1935) with Bakun & Joyner (1984) for MEL and ADE magnitudes post 2002 and 2004, respectively and correct to MLM92
	- use MEL corrections for all RC solutions
	- assume Ricther (1935) for all QEDB solutions and correct to MLM92
	- use Russell Cuthbertson's updated list of Australian seismic stations to determine event-station distances
3) Run Preferred\reviseMLs_2017.m:
	- generates mdat_pref12.mat
4) Run "Append_Preferred_MWs\refmt_au_mw_tables.py"
	- makes "combined_au_mw.dat"
5) Run "Append_Preferred_MWs\append_mw.m"
	- append MW estimates from "combined_au_mw.dat"
	- generates mdat_mw_pref12.mat
6) Run "Preferred\get_pref_MW.m"
	- make final catalogue
