# Minimal example
python shapefile_2_complexfault.py -shapefile SUMATRA/sumatra.shp -shapefile_depth_attribute level

# More complex
# Complex fault example
python shapefile_2_complexfault.py -shapefile SUMATRA/sumatra.shp -shapefile_depth_attribute level -source_model_name SUMATRA_1234 -complex_fault_id 1234 -complex_fault_name sumatra_GAR_version -complex_fault_tectonic_region SUMATRA -magnitude_scaling_relation WC1994 -rupture_aspect_ratio 2.0 -a_value 20 -b_value 1.0 -min_mag 4.5 -max_mag 9.7 -rake 90 -output_dir nrml

# Simple fault example
python shapefile_2_simplefault.py -shapefile ADELAIDE/Adelaide_simple_faults.shp -shapefile_faultname_attribute Name -shapefile_dip_attribute Dip -shapefile_sliprate_attribute Slip_Rate -source_model_name Adelaide_faults -simple_fault_tectonic_region Non_cratonic -magnitude_scaling_relation WC1994 -rupture_aspect_ratio 1.5 -upper_depth 0.001 -lower_depth 20 -a_value 1.0 -b_value 1.0 -min_mag 4.5 -max_mag 7.0 -rake 90 -output_dir Adelaide_test

# Simple fault CE example (note max_mag here means characteristic magnitude)
python shapefile_2_simplefault_CE.py -shapefile ADELAIDE/Adelaide_simple_faults.shp -shapefile_faultname_attribute Name -shapefile_dip_attribute Dip -shapefile_sliprate_attribute Slip_Rate -source_model_name Adelaide_faults -simple_fault_tectonic_region Non_cratonic -magnitude_scaling_relation WC1994 -rupture_aspect_ratio 1.5 -upper_depth 0.001 -lower_depth 20 -a_value 1.0 -b_value 1.0 -min_mag 4.5 -max_mag 7.0 -rake 90 -output_dir Adelaide_test_CE -incremental_mfd True

# Simple fault MM example (maximum magnitude model)
python shapefile_2_simplefault_MM.py -shapefile faults/cadell_fault.shp -shapefile_faultname_attribute Name -shapefile_dip_attribute Dip -shapefile_sliprate_attribute SL_RT_LT -source_model_name Cadell_fault -simple_fault_tectonic_region Non_cratonic -magnitude_scaling_relation WC1994 -rupture_aspect_ratio 1.5 -upper_depth 0.001 -lower_depth 20 -a_value 1.0 -b_value 1.0 -min_mag 4.5 -max_mag 7.0 -rake 90 -output_dir Cadell_MM -incremental_mfd True
