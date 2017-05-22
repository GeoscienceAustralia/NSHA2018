# coding: utf-8

#get_ipython().magic(u'matplotlib inline')

from hmtk.seismicity.smoothing.smoothed_seismicity import SmoothedSeismicity
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser

# Python dependences
import os, sys
import h5py
import numpy as np   # Numpy - Python's numerical library
import matplotlib.pyplot as plt  # Matplotlib - Python's plotting library
from copy import deepcopy   # Python module for copying objects

# Input and Output Tools
# Catalogue and sources 
from hmtk.parsers.catalogue import CsvCatalogueParser   # Reads an earthquake catalogue from CSV
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueWriter  # Writes an earthquake catalogue to CSV
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser  # Imports a source model from XML

# Plotting tools
#from hmtk.plotting.mapping import HMTKBaseMap
from hmtk.plotting.seismicity.completeness import plot_stepp_1972
from hmtk.plotting.seismicity.catalogue_plots import plot_magnitude_time_scatter
from hmtk.plotting.seismicity.catalogue_plots import plot_depth_histogram
from hmtk.plotting.seismicity.catalogue_plots import plot_magnitude_time_density
from hmtk.plotting.seismicity.max_magnitude.cumulative_moment import plot_cumulative_moment 
from hmtk.plotting.seismicity.catalogue_plots import (plot_observed_recurrence, 
                                                      get_completeness_adjusted_table,
                                                      _get_catalogue_bin_limits)

# Seismicity tools: Events and declustering methods
from hmtk.seismicity.selector import CatalogueSelector
from hmtk.seismicity.declusterer.dec_afteran import Afteran 
from hmtk.seismicity.declusterer.dec_gardner_knopoff import GardnerKnopoffType1 
from hmtk.seismicity.declusterer.distance_time_windows import (GardnerKnopoffWindow, 
                                                               GruenthalWindow, 
                                                               UhrhammerWindow)

# Completeness tools
from hmtk.seismicity.completeness.comp_stepp_1971 import Stepp1971

# Seismicity tools: Recurrence methods
from hmtk.seismicity.occurrence.aki_maximum_likelihood import AkiMaxLikelihood
from hmtk.seismicity.occurrence.b_maximum_likelihood import BMaxLikelihood
from hmtk.seismicity.occurrence.kijko_smit import KijkoSmit
from hmtk.seismicity.occurrence.weichert import Weichert

# Seismicity tools: Recurrence methods
from hmtk.seismicity.max_magnitude.kijko_sellevol_fixed_b import KijkoSellevolFixedb
from hmtk.seismicity.max_magnitude.kijko_sellevol_bayes import KijkoSellevolBayes
from hmtk.seismicity.max_magnitude.kijko_nonparametric_gaussian import KijkoNonParametricGaussian
from hmtk.seismicity.max_magnitude.cumulative_moment_release import CumulativeMoment 

# Seismicity tools: Smoothed seismicity
from hmtk.seismicity.smoothing.smoothed_seismicity import SmoothedSeismicity 
from hmtk.seismicity.smoothing.kernels.isotropic_gaussian import IsotropicGaussian 

#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueWriter
import helmstetter_werner_2012 as h_w
# To build source model
from hmtk.sources.source_model import mtkSourceModel
from hmtk.sources.point_source import mtkPointSource
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.pmf import PMF
print "Everything Imported OK!"

bvalue = float(sys.argv[1])
print 'b value', bvalue

ifile = "../../catalogue/data/AUSTCAT_V0.12_hmtk_declustered.csv"
#ifile = "../../catalogue/data/AUSTCAT_V0.12_hmtk_mx_orig.csv"
parser = CsvCatalogueParser(ifile)
catalogue = parser.read_file(start_year=1965, end_year=2010)
# How many events in the catalogue?
print "The catalogue contains %g events" % catalogue.get_number_events()
neq = len(catalogue.data['magnitude'])
print "The catalogue contains %g events" % neq
# What is the geographical extent of the catalogue?
bbox = catalogue.get_bounding_box()
print "Catalogue ranges from %.4f E to %.4f E Longitude and %.4f N to %.4f N Latitude\n" % bbox

catalogue.sort_catalogue_chronologically()
index = np.logical_and(catalogue.data["magnitude"] > 1.5, catalogue.data["depth"] >= 0.0) 
catalogue.purge_catalogue(index)
catalogue.get_number_events()

# Copying the catalogue and saving it under a new name "catalogue_clean"
catalogue_clean = deepcopy(catalogue)
# remove nan magnitudes
catalogue_clean.sort_catalogue_chronologically()
catalogue_clean.data['magnitude']
catalogue_clean.data['year']
catalogue_clean.get_decimal_time()
catalogue_clean.data['longitude']

catalogue_depth_clean = deepcopy(catalogue_clean)
index = catalogue_depth_clean.data['depth']>=0.
catalogue_depth_clean.purge_catalogue(index)
catalogue_clean.get_number_events()

#source_model_file = "../zones/2012_mw_ge_4.0/NSHA13_Background/input/best/NSHA13_BACKGROUND_best.xml"
source_model_file = 'Aus_cont_testzone.xml'
parser = nrmlSourceModelParser(source_model_file)

# Parse the seismic sources and save them into a variable called "source_model"
source_model = parser.read_file("Aus Source Model 1") # You need to supply a name for the source model

# Map configuration
llon, ulon, llat, ulat = catalogue_clean.get_bounding_box()
#map_config = {'min_lon': np.floor(llon), 'max_lon': np.ceil(ulon),
 #             'min_lat': np.floor(llat), 'max_lat': np.ceil(ulat), 'resolution':'c'}
map_config = {'min_lon': np.floor(100), 'max_lon': np.ceil(160),
              'min_lat': np.floor(-45), 'max_lat': np.ceil(-4), 'resolution':'c'}
# Creating a basemap - input a cconfiguration and (if desired) a title
#basemap1 = HMTKBaseMap(map_config, 'Earthquake Catalogue')

# Adding the seismic sources
#basemap1.add_source_model(source_model, area_border='r-', border_width=1.5, alpha=0.5)

# Select catalogue from within sourcezone
selector1 = CatalogueSelector(catalogue_depth_clean, create_copy=True)
for source in source_model.sources:
    source.select_catalogue(selector1)
    
    llon, ulon, llat, ulat = source.catalogue.get_bounding_box()
    print llon, ulon, llat, ulat
    # Map the Source
#    src_basemap = HMTKBaseMap(map_config, "Source: {:s}".format(source.name))
 #   print "Source ID: %s  Source Name: %s   Number of Events: %g" % (source.id, source.name,
 #                                                                    source.catalogue.get_number_events())
    # Add on the catalogue
 #   src_basemap.add_catalogue(source.catalogue, overlay=False)
#completeness_table_a = np.array([[1965., 4.0]])
#completeness_table_a = np.array([[1980., 3.5],
#                                 [1965., 4.0]])
completeness_table_a = np.array([[1990., 3.0],
                                 [1980., 3.5],
                                 [1965., 4.0]])
###plot_magnitude_time_density(source.catalogue, 0.1, 1.0,
#                            completeness=completeness_table_a)

grid_lims = [110., 160.0, 0.1, -45.0, -5.0, 0.1, 0., 20., 20.]

try:
    os.remove("Aus1_tmp.hdf5")
except OSError:
    pass
config = {"bandwidth": 50,
          "r_min": 1.0E-7, 
          "bvalue": bvalue, "mmin": 3.0,
          "learning_start": 1965, "learning_end": 2003,
          "target_start": 2004, "target_end": 2010}
#smoother = h_w.FixedWidthSmoothing(grid_lims, config, source.catalogue, storage_file="Aus1_tmp.hdf5")
#smoother._get_catalogue_completeness_weights(completeness_table_a)
#smoother.build_catalogue_2_grid_array()
#smoother.run_smoothing(config["r_min"], config["bandwidth"])

#np.savetxt("Australia_FixedWidth2_50km.csv",
 #          np.column_stack([smoother.grid, smoother.rates]),
 #          delimiter=",",
 #          fmt=["%.4f", "%.4f", "%.8e"],
 #         header="longitude,latitude,rate" 
 #         )

try:
    os.remove(("Aus1_tmp2%.3f.hdf5" % bvalue))
except OSError:
    pass
config = {"k": 4,
          "r_min": 1.0E-6, 
          "bvalue": bvalue, "mmin": 3.0,
          "learning_start": 1965, "learning_end": 2009,
          "target_start": 2009, "target_end": 2010}

smoother = h_w.HelmstetterEtAl2007(grid_lims, config, source.catalogue, storage_file=("Aus1_tmp2%.3f.hdf5" % bvalue))

smoother._get_catalogue_completeness_weights(completeness_table_a)
smoother.build_distance_arrays()
smoother.build_catalogue_2_grid_array()
# Exhaustive smoothing
exhaustive = False
if exhaustive == True:
    params, poiss_llh = smoother.exhaustive_smoothing(np.arange(3,6,1), np.arange(1.0e-6,1.0e-5,2.0e-6))
    print params, poiss_llh
    smoother.config["k"] = params[0]
    smoother.config["r_min"] = params[1]
d_i = smoother.optimise_bandwidths()
smoother.run_smoothing(config["r_min"], d_i)
smoother_filename = "Australia_Adaptive_K%i_b%.3f_mmin%.1f.csv" % (smoother.config['k'], smoother.config['bvalue'], smoother.config['mmin'])
np.savetxt(smoother_filename,
           np.column_stack([smoother.grid, smoother.rates]),
           delimiter=",",
           fmt=["%.4f", "%.4f", "%.8e"],
           header="longitude,latitude,rate" 
           )

# Creating a basemap - input a cconfiguration and (if desired) a title
#basemap1 = HMTKBaseMap(map_config, 'Smoothed seismicity rate')
#basemap1.m.drawmeridians(np.arange(llat, ulat, 5))
#basemap1.m.drawparallels(np.arange(llon, ulon, 5))
# Adding the smoothed grip to the basemap
#sym = (2., 3.,'cx')
#x,y = basemap1.m(smoother.grid[:,0], smoother.grid[:,1])
#basemap1.m.scatter(x, y, marker = 's', c = np.log10(smoother.rates), cmap = plt.cm.coolwarm, zorder=10, lw=0, vmin=-6.5, vmax=1.5)
#basemap1.m.drawcoastlines(linewidth=1, zorder=50) # Add coastline on top
#basemap1.m.drawmeridians(np.arange(llat, ulat, 5))
#basemap1.m.drawparallels(np.arange(llon, ulon, 5))
#plt.colorbar()#label='log10(Smoothed rate per cell)')
#plt.legend()
#basemap1.m.scatter(x, y, marker = 's', c = smoother.data[:,4], cmap = plt.cm.coolwarm, zorder=10)
#basemap1.m.scatter([150],[22], marker='o')
#basemap1.fig.show()

#(smoother.data[0], smoother.data[1])
#basemap1.add_catalogue(catalogue_depth_clean, erlay=False)
#plt.savefig('./smoothed_%i_%i_%i_%i_%i_%.6f.png' %             (smoother.config["k"], smoother.config["learning_start"],
#            smoother.config["learning_end"], smoother.config["target_start"],
#            smoother.config["target_end"],smoother.config["r_min"]))

# Build nrml input file of point sources

source_list = []
#i=0
min_mag = 4.5
max_mag = 7.2
#bval = 1.0 # just define as 1 for time being
# Read in data again to solve number fomatting issue in smoother.data
# For some reason it just returns 0 for all a values
data = np.genfromtxt(smoother_filename, delimiter = ',', skip_header = 1)
#print max(data[:,4])
#print data[:,4]
#print len(data[:,4])
tom = PoissonTOM(50) # Dummy temporal occurence model for building pt sources
msr = Leonard2014_SCR()
for j in range(len(data[:,2])):
#    print smoother.data[j,:]
    identifier = 'ASS' + str(j)
    name = 'Helmstetter' + str(j)
    point = Point(data[j,0],data[j,1],
                10)
    rate = data[j,2]
    # Convert rate to a value??
    #aval = rate
    aval = np.log10(rate) + bvalue*config["mmin"]
    #print rate, aval
    mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, aval, bvalue)
    hypo_depth_dist = PMF([(0.5, 10.0),
                          (0.25, 5.0),
                          (0.25, 15.0)])
    nodal_plane_dist = PMF([(0.3, NodalPlane(0, 30, 90)),
                            (0.2, NodalPlane(90, 30, 90)),
                            (0.3, NodalPlane(180, 30, 90)),
                            (0.2, NodalPlane(270, 30, 90))])
    point_source = PointSource(identifier, name, 'Non_cratonic',
                               mfd, 2, msr,
                               2.0, tom, 0.1, 20.0, point,
                               nodal_plane_dist, hypo_depth_dist)
#                           upper_depth = 0.1, lower_depth = 20.0,
#                           nodal_plane_dist = nodal_plane_dist,
#                           hypo_depth_dist = hypo_depth_dist)
    source_list.append(point_source)
#    i+=1
#    if j==1000:
#        break
filename = "Australia_Adaptive_K%i_b%.3f_mmin%.1f.xml" % (smoother.config['k'], smoother.config['bvalue'], smoother.config['mmin'])
mod_name = "Australia_Adaptive_K%i_b%.3f" % (smoother.config['k'], smoother.config['bvalue'])   
nodes = list(map(obj_to_node, sorted(source_list)))
source_model = Node("sourceModel", {"name": name}, nodes=nodes)
with open(filename, 'wb') as f:
    nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
#source_model = mtkSourceModel(identifier=0, name=mod_name,
#                              sources = source_list)
#source_model.serialise_to_nrml('source_model_' + mod_name + '.xml')

