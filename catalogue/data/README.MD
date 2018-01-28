List of catalogues used to determine MFDs for the 2018 NSHA:
	
* NSHA18CAT.MW.V0.1.csv: Raw Australian earthquake catalogue developed for the 2018 NSHA:
	- Preserves all magnitude types and converted magnitudes
	- Uses hierarchical structure as outlined in Allen et al (in prep) to determine preferred magnitude type
	- Uses Ghasemi & Allen (2017) as preferred ML-MW conversions
	- Uses updated mb-MW and MS-MW conversions as outlined in Allen et al (in prep)
	- Preferred magnitude column is MX_REVML
	- Further details on the catalogue's development is [here](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/catalogue/data/2012_catalogue_readme.txt) (to be updated)

* NSHA18CAT_V0.1_hmtk_mx_orig.csv: full catalogue developed for the 2018 NSHA in [Hazard Modellers Toolkit (HMTK)](https://github.com/GEMScienceTools/hmtk) format:
	- Based on NSHA18CAT.MW.V0.1.csv
	- Preferred magnitude is [MX_REVML](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt) as required for the Leonard (2008) declustering algorithm
	- Catalogue adds several additional magnitude columns to standard HMTK format.  Details are [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt) and preferred Mw as per 2018 conversions
	- Not declustered
	
* NSHA18CAT_V0.1_hmtk_declustered.csv: declustered catalogue developed for the 2018 NSHA in HMTK format:
	- Based on NSHA18CAT.MW.V0.1.csv
	- Preferred magnitude is in MW as per 2018 conversions
	- Declustered using Leonard (2008) algorithm
	- Catalogue adds several additional magnitude columns to standard HMTK format.
	
* Catalogues used to generate magnitude conversion equations (Allen et al, in prep):
	- NSHA18CAT.ML-MW.csv
	- NSHA18CAT.mb-MW.csv
	- NSHA18CAT.MS-MW.csv