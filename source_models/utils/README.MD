utils
=============

Code for manipulating source models files:
- Covnerting formats:
  - area2pt_source.py converts an area source into its component point sources (as OpenQuake would do)
  - nrml2shp.py converts OpenQuake nrml to a shapefile. 
  `python nrml2shp.py <source_model_nrml_file.xml>`
- Spatial filering
  - pt2fault_distance.py identifies ruptures in a point source model that intersect a fault source model and elimiates these
