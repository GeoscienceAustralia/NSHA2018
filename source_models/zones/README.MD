NSHA18 source model workflow
============================

*	Base shapefiles are developed [here](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/source_models/zones/shapefiles) from original data sources:
	- Bespoke scripts are written to convert each area source model into a shapefile with common attributes, [for example](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/shapefiles/ARUP/make_ARUP_model.py)
	- Where possible, attributes are filled, otherwise null values are used
	- Base shapefiles generated have the suffix *_NSHA18.shp

* The script [get_NSHA18_mfds_merge.py](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/get_NSHA18_mfds_merge.py)  uses these files to generate earthquake magnitude-frequency statistics:
	- The script parses a parameter file of [this](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/2018_mw/ARUP/arup.param) format
	- The parameter file provides a path to the base shapefiles above
	- For each area source, magnitude-frequency statistics are calculated based on the catalogue specified in the shapefile
	- The script outputs:
		- Zone specific statistics, [for example](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/source_models/zones/2018_mw/ARUP/mfd)
		- [Summary files](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/source_models/zones/2018_mw/ARUP), appended with [*_NSHA18.MERGE.pdf](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/2018_mw/ARUP/ARUP_NSHA18_MFD.MERGE.pdf) and [*_NSHA18_MFD.csv](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/2018_mw/ARUP/ARUP_NSHA18_MFD.csv)
		- Output [shapefiles](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/source_models/zones/2018_mw/ARUP/shapefiles) with completed magnitude-frequency statistics
		
* The script [write_oq_inputs.py](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/zones/write_oq_inputs.py) parses these newly processed shapefiles to output OpenQuake-engine [input files](https://github.com/GeoscienceAustralia/NSHA2018/tree/master/source_models/zones/2018_mw/AUS6/input)

* The script [make_source_model_logic_tree.py](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/logic_trees/make_source_model_logic_tree.py) builds the OpenQuake-engine [source model logic tree](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/source_models/logic_trees/nsha18_source_model_logic_tree.xml) and copies the necessary source model files to the appropriate directory for running the combined hazard model:
	- Source model weights are determined from the [results](https://github.com/GeoscienceAustralia/NSHA2018/blob/master/shared/seismic_source_model_weights_rounded_p0.4.edit.csv) of the Seismic Source Characterisation expert elicitation workshop
	