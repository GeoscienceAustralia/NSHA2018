List of catalogues used to determine MFDs for the 2018 NSHA:
	
* NSHA18CAT.MW.V0.1.csv: Raw Australian earthquake catalogue developed for the 2018 NSHA:
	- Preserves all magnitude types and converted magnitudes
	- Uses hierarchical structure as outlined in Allen et al (in prep) to determine preferred magnitude type
	- Uses Ghasemi & Allen (2017) as preferred ML-MW conversions
	- Uses updated mb-MW and MS-MW conversions as outlined in Allen et al (in prep)
	- Preferred magnitude column is MX_REVML
	- Further details on the catalogue's development is [here](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/info/NSHA18_description_of_files_and_metadata.XLS) (to be updated)

* NSHA18CAT_V0.1_hmtk_mx_orig.csv: full catalogue developed for the 2018 NSHA in [Hazard Modellers Toolkit (HMTK)](https://github.com/GEMScienceTools/hmtk) format:
	- Based on NSHA18CAT.MW.V0.1.csv
	- Preferred magnitude is [mx_origML](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt) as required for the Leonard (2008) declustering algorithm
	- Catalogue adds several additional magnitude columns to standard HMTK format.  Details are [here](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/info/NSHA18_description_of_files_and_metadata.XLS) and preferred Mw as per 2018 conversions
	- Not declustered
	
* NSHA18CAT_V0.1_hmtk_declustered.csv: declustered catalogue developed for the 2018 NSHA in HMTK format:
	- Based on NSHA18CAT.MW.V0.1.csv
	- Declustered using Leonard ([2008](https://pubs.geoscienceworld.org/ssa/bssa/article/98/3/1458/341978)) algorithm using [mx_origML](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/info/NSHA18_description_of_files_and_metadata.XLS)
	- Preferred magnitude is in MW as per 2018 conversions
	- Catalogue adds several additional magnitude columns to standard HMTK format.
	
* Catalogues used to generate magnitude conversion equations (Allen et al, in prep):
	- NSHA18CAT.ML-MW.csv
	- NSHA18CAT.mb-MW.csv
	- NSHA18CAT.MS-MW.csv
	
* ISC-GEM_V5_hmtk_full.csv: full [International Seismological Centre-Global Earthquake Model (ISC-GEM) earthquake catalogue](http://www.isc.ac.uk/iscgem/overview.php) in HMTK format:
	- Catalogue to be used for adjacent offshore sources
	- Refomated from original catalogue format as described in Storchak et al. ([2015](https://www.sciencedirect.com/science/article/pii/S003192011400154X))
	- Version 5 released on 2018-02-27 for the period of 1904-2014
	
* ISC-GEM_V5_hmtk_GK74_declustered.csv: declustered ISC-GEM earthquake catalogue in HMTK format:
	- Catalogue as described above, decustered using the Gardner and Knopoff ([1974](https://pubs.geoscienceworld.org/ssa/bssa/article/64/5/1363/117341/is-the-sequence-of-earthquakes-in-southern)) algorithm

* ISC-GEM_V4_hmtk_full.csv: full International Seismological Centre-Global Earthquake Model (ISC-GEM) earthquake catalogue in HMTK format:
	- Catalogue to be used for adjacent offshore sources
	- Refomated from original catalogue format as described in Storchak et al. ([2015](https://www.sciencedirect.com/science/article/pii/S003192011400154X))
	- Version 4 released on 2017-01-26 for the period of 1900-2013
	
* ISC-GEM_V4_hmtk_GK74_declustered.csv: declustered ISC-GEM earthquake catalogue in HMTK format:
	- Catalogue as described above, decustered using the Gardner and Knopoff ([1974](https://pubs.geoscienceworld.org/ssa/bssa/article/64/5/1363/117341/is-the-sequence-of-earthquakes-in-southern)) algorithm
	
* GIS files for catalogues are provided [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/shapefiles)