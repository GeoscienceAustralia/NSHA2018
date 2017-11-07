utils
=============

Code for manipulating source models files:
- Converting formats:
  - area2pt_source.py converts an area source into its component point sources (as OpenQuake would do)
  - nrml2shp.py converts OpenQuake nrml to a shapefile. 
  `python nrml2shp.py <source_model_nrml_file.xml>`
- Spatial filtering
  - pt2fault_distance.py identifies ruptures in a point source model that intersect a fault source model and elimiates these
  - pt_in_polygon.py identifies point sources that are within a polygon (defined at present by an area source zone) and performs some action (e.g. modifying a parameter based on which zone the point is in)
