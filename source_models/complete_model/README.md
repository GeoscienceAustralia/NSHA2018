Steps to build NSHA18 source models on the NCI
==============================================    

1) run "write_oq_inputs.py" for all area sources:
	- this makes OQ input files and adds Java/Banda Sea/PNG fault-source model to all area sources
	- File appended is "../banda/Banda_Fault_Sources_NSHA_2018.xml"
	
2) Run "append_banda_nfsm.py" for all smoothed seismicity sources:
	- this takes the default smoothed seismicity models and appends Java/Banda Sea/PNG area and fault-source models
	- File appended is "../zones/2018_mw/Java_Banda_PNG/input/collapsed/Java_Banda_PNG_collapsed.xml"
	- script optionally appends the national fault source model
	- File appended is "../faults/NFSM/National_Fault_Source_Model_2018_Collapsed_all_methods_collapsed_inc_cluster.xml"
	
3) Untar seismotectonic models

4) Run "make_source_model_logic_tree.py"
	- copies all input source models from disparate directories to final run directory
	- builds OQ source model logic tree (needs to be copied to working directory)