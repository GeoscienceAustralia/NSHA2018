List of catalogues used to determine MFDs for the 2018 NSHA:
	
* AUSTCAT.MW.V0.11.csv: Raw Australian earthquake catalogue developed for the 2012 NSHA:
	- Preserves all magnitude types and converted magnitudes
	- Uses hierarchical structure as outlined in Appendix A of Burbidge (2012) to determine preferred magnitude type
	- Uses Grünthal et al. (2003) ML-MW conversions
	- Preferred magnitude column is MX_REVML
	- Further details on the catalogue's development is [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt)

* AUSTCAT_V0.11_hmtk.csv: full catalogue developed for the 2012 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.11.csv
	- Preferred magnitude is in MW as per 2012 conversions
	- Not declustered
	
* AUSTCAT_V0.11_hmtk_declustered.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.11.csv
	- Preferred magnitude is in MW as per 2012 conversions
	- Declustered using Leonard (2008) algorithm
	
* AUSTCAT_V0.11_hmtk_declustered_GK74.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Preferred magnitude is in MW as per 2012 conversions
	- Declustered using Gardner & Knopoff (1974) algorithm, as implemented in the HMTK

* AUSTCAT_V0.11_hmtk_mx_revml.csv: full catalogue developed for the 2012 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.11.csv
	- Uses hierarchical structure as outlined in Appendix A of Burbidge (2012) to determine preferred magnitude type
	- Preferred magnitude is [MX_REVML](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt)
	- Not declustered
	
* AUSTCAT_V0.11_hmtk_mx_revml_declustered.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.11.csv
	- Uses hierarchical structure as outlined in Appendix A of Burbidge (2012) to determine preferred magnitude type
	- Preferred magnitude is [MX_REVML](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt)
	- Declustered using Leonard (2008) algorithm

* AUSTCAT.MW.V0.12.csv: Raw Australian earthquake catalogue updated for initial use for the 2018 NSHA:
	- Based on AUSTCAT.MW.V0.11.csv
	- Updated logic for ML corrections
	- Includes Ghasemi et al (2016) MW estimates where applicable
	- Uses updated 2017 ML-MW conversions based on existing Australian data
	- Uses GEM-ISC MS-MW conversions (Di Giacomo et al., 2015) where applicable
	- Further details [here](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/catalogue/matlab/README.MD)
	- Background on the 2012 catalogue's development and a description of fields is [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt)
	
* AUSTCAT_V0.12_hmtk.csv: full catalogue developed for the 2018 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.12.csv
	- Preferred magnitude is in MW as per 2017 conversions
	- Not declustered
	
* AUSTCAT_V0.12_hmtk_declustered.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Based on AUSTCAT.MW.V0.12.csv
	- Preferred magnitude is in MW as per 2017 conversions
	- Declustered using Leonard (2008) algorithm
