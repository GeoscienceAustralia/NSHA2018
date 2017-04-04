List of catalogues used to determine MFDs for the 2018 NSHA:

* AUSTCAT_hmtk.csv: full catalogue developed for the 2012 NSHA in HMTK format:
	- Preferred magnitude is in MW as per 2012 conversions
	- Not declustered
	
* AUSTCAT_hmtk_declustered.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Preferred magnitude is in MW as per 2012 conversions
	- Declustered using Leonard (2008) algorithm
	
* AUSTCAT_hmtk_declustered_GK74.csv: declustered catalogue developed for the 2012 NSHA in HMTK format:
	- Preferred magnitude is in MW as per 2012 conversions
	- Declustered using Gardner & Knopoff (1974) algorithm, as implemented in the HMTK
	
* AUSTCAT.MW.V0.11.csv: Raw Australian earthquake catalogue developed for the 2012 NSHA:
	- Preserves all magnitude types and converted magnitudes
	- Uses hierarchical structure as outlined in Appendix A of Burbidge (2012) to determine preferred magnitude type
	- Uses Grünthal et al. (2003) ML-MW conversions
	- Preferred magnitude column is MX_REVML
	- Further details on the catalogue's development is [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/catalogue/data/2012_catalogue_readme.txt)