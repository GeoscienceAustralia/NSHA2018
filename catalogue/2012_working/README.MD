The 2012 NSHA catalogue is used as the starting poitnt for the development of the draft design hazard values as prepared for the 2019 AS1170.4.  Geoscience Australia intend to update the catalogue prior to the finalisation of seismic design values in late 2017.  The steps used to develop the catalogue are:

1) Use 2012 Catalogue V0.11 as the starting point (see Burbidge, 2012 Appendix A)
2) Update ML correction code to:
	- replace Richter (1935) with Bakun & Joyner (1984) for MEL and ADE magnitudes post 2002 and 2004, respectively and correct to MLM92
	- use MEL corrections for all RC solutions
	- assume Ricther (1935) for all QEDB solutions and correct to MLM92
	- use Russell Cuthbertson's updated list of Australian seismic stations to determine event-station distances
3) Run "matlab\peferred\reviseMLs_2017.m":
	- generates "mdat_pref12.mat"
4) Run "matlab\append_mw\refmt_au_mw_tables.py"
	- makes "combined_au_mw.dat"
5) Run "append_mw\append_mw.m"
	- append MW estimates from "combined_au_mw.dat"
	- generates "mdat_mw_pref12.mat"
6) Run "matlab\preferred\get_pref_MW.m"
	- make final catalogue: "data\AUSTCAT.MW.V0.12.csv"
7) A number of duplicate events were identified in the catalogue.  These were removed by running: "remove_duplicate_events.py"
9) Declustering was performed using "decluster.py":
	- Deculsters using Leonard (2008) algorithm on original catalogue magnitudes (e.g., ML, MS, mb, etc.).  This declustering algorithm was developed using original catalogue magnitudes, so this is appropriate.
	- Deculstering performed on "AUSTCAT_V0.12_hmtk_mx_orig.csv"
	- Generates "AUSTCAT_V0.12_hmtk_declustered.csv" and supplants preferred MW into "magnitude" column after declustering
