Contains files for building the Australian smoothed seismicity models.
Jobs 1 and 2 can be run concurrently; once they have finished jobs 3 and 4 can also be run concurrently.

1. adaptive_smoothing_parallel.py (run with run_adaptive_parallel.sh on NCI) 
This script will extract combinations of bvalues and completeness models from the 'domains' shapefile specified by the user defined parameter domains_shp. It is currently using pre-optimsed values of K and r_min or the Helmstetter et al 2007 method (i.e. these values are not separately optimsed for each bvalue-completeness combination). 

2. fixed_smoothing_parallel.py (run with run_frankel_parallel.sh on NCI). 
Implements the frankel smoothing method for each combination of bvalues and completeness models from the 'domains' shapefile specified by the user defined parameter domains_shp.

3. combine_ss_models.py (run with run_combine_ss.sh on NCI)
Uses the domains shapefiles to filter points from the relevant smoothed seismicity models. It then merges rates for the different bvalues (best, upper, lower) and their given weights. 

4. combine_ss_frankel_models.py (run with run_combine_ss_frankel.sh on NCI) 
As above, but pointing to the relevant fixed smoothing models.