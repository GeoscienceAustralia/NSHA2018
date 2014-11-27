# Minimal example
python shapefile_2_complexfault.py -shapefile SUMATRA/sumatra.shp -shapefile_depth_attribute level

# More complex
python shapefile_2_complexfault.py -shapefile SUMATRA/sumatra.shp -shapefile_depth_attribute level -source_model_name SUMATRA_1234 -complex_fault_id 1234 -complex_fault_name sumatra_GAR_version -complex_fault_tectonic_region SUMATRA -magnitude_scaling_relation WC1994 -rupture_aspect_ratio 2.0 -a_value 20 -b_value 1.0 -min_mag 4.5 -max_mag 9.7 -rake 90
