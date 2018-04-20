# coding: utf-8
from hmtk.seismicity.smoothing.smoothed_seismicity import SmoothedSeismicity

# Python dependences
import os, sys
import h5py
import numpy as np   # Numpy - Python's numerical library
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # Matplotlib - Python's plotting library
from copy import deepcopy   # Python module for copying objects
import ogr
import shapefile
from shapely.geometry import Point, Polygon
from utilities import params_from_shp

# For running in parallel
import time
from time import localtime, strftime, gmtime
import string
import pypar

# Input and Output Tools
# Catalogue and sources 
from hmtk.parsers.catalogue import CsvCatalogueParser   # Reads an earthquake catalogue from CSV
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueWriter  # Writes an earthquake catalogue to CSV
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser  # Imports a source model from XML

# Plotting tools
from hmtk.plotting.mapping import HMTKBaseMap
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

#bvalue = float(sys.argv[1])
#print 'b value', bvalue
#domains_shp = '../zones/2018_mw/Domains_mc_ext/shapefiles/Domains_NSHA18_MFD.shp'
domains_shp = '../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
ifile = "../../catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv"
#ifile = "../../catalogue/data/AUSTCAT_V0.12_hmtk_mx_orig.csv"
# Flag for whether to overwrite exiting .xml source model 
# files with the same b value and completeness combination.
# Shoudld normally set to True unless you are being really careful.
overwrite = True
#####################################
# Shouldn't need input below here
#####################################

def run_smoothing(grid_lims, config, catalogue, completeness_table,map_config, run):
    """Run all the smoothing
    :params config:
        Dictionary of configuration parameters.
        For more info see helmstetter_werner_2012 code 
        and docs.
    """

    completeness_string = 'comp'
    for ym in completeness_table:
        completeness_string += '_%i_%.1f' % (ym[0], ym[1])
    smoother_filename = "Australia_Adaptive_K%i_b%.3f_mmin%.1f_%s.csv" % (
        config['k'], config['bvalue'], config['mmin'],
        completeness_string)

    filename = smoother_filename[:-4] + '.xml'
    if os.path.exists(filename) and not overwrite:
        print '%s already created, not overwriting!' % filename
        return

    smoother = h_w.HelmstetterEtAl2007(grid_lims, config, catalogue, 
                                       storage_file=("Aus1_tmp2%.3f_%s.hdf5" % (config['bvalue'],run)))
    smoother._get_catalogue_completeness_weights(completeness_table)
    smoother.build_distance_arrays()
    smoother.build_catalogue_2_grid_array()
    # Exhaustive smoothing
    exhaustive = False
    if exhaustive == True:
        params, poiss_llh = smoother.exhaustive_smoothing(np.arange(2,10,1), np.arange(1.0e-6,1.0e-5,2.0e-6))
        print params, poiss_llh
        smoother.config["k"] = params[0]
        smoother.config["r_min"] = params[1]
    #print 'Exiting now, re-run using optimised parameters'
    #sys.exit()
    d_i = smoother.optimise_bandwidths()
    smoother.run_smoothing(config["r_min"], d_i)
    data =  np.column_stack([smoother.grid, smoother.rates])
    np.savetxt(smoother_filename,
               data,
#               np.column_stack([smoother.grid, smoother.rates]),
               delimiter=",",
               fmt=["%.4f", "%.4f", "%.8e"],
               header="longitude,latitude,rate" 
               )
                                       
    # Creating a basemap - input a cconfiguration and (if desired) a title
    title = 'Smoothed seismicity rate for learning \nperiod %i %i, K=%i, Mmin=%.1f' % (
        config['learning_start'], config['learning_end'], smoother.config['k'], smoother.config['mmin'])
    basemap1 = HMTKBaseMap(map_config, title)
    basemap1.m.drawmeridians(np.arange(map_config['min_lat'],
                                       map_config['max_lat'], 5))
    basemap1.m.drawparallels(np.arange(map_config['min_lon'],
                                        map_config['max_lon'], 5))
    # Adding the smoothed grip to the basemap
    sym = (2., 3.,'cx')
    x,y = basemap1.m(smoother.grid[:,0], smoother.grid[:,1])
    if smoother.config['mmin'] == 3.5:
        vmax=-1.0
    elif smoother.config['mmin'] == 4.0:
        vmax=-2.5
    else:
        vmax=-1.0
    basemap1.m.scatter(x, y, marker = 's', c = np.log10(smoother.rates), cmap = plt.cm.coolwarm, zorder=10, lw=0, vmin=-7.0, vmax=vmax)
    basemap1.m.drawcoastlines(linewidth=1, zorder=50) # Add coastline on top
    #basemap1.m.drawmeridians(np.arange(llat, ulat, 5))
    #basemap1.m.drawparallels(np.arange(llon, ulon, 5))
    plt.colorbar(label='Log10(Smoothed rate per cell)')
    #plt.colorbar()#label='log10(Smoothed rate per cell)')
    plt.legend()
    #basemap1.m.scatter(x, y, marker = 's', c = smoother.data[:,4], cmap = plt.cm.coolwarm, zorder=10)
    #basemap1.m.scatter([150],[22], marker='o')
    #basemap1.fig.show()

    #(smoother.data[0], smoother.data[1])
    #basemap1.add_catalogue(catalogue_depth_clean, erlay=False)
    figname = smoother_filename[:-4] + '_smoothed_rates_map.png'
    plt.savefig(figname)
                                       
    source_list = []
    #i=0
    min_mag = 4.5
    max_mag = 7.2
    # Read in data again to solve number fomatting issue in smoother.data
    # For some reason it just returns 0 for all a values
    #data = np.genfromtxt(smoother_filename, delimiter = ',', skip_header = 1)

    tom = PoissonTOM(50) # Dummy temporal occurence model for building pt sources
    msr = Leonard2014_SCR()
    for j in range(len(data[:,2])):
        identifier = 'ASS' + str(j) + '_' + str(run)
        name = 'Helmstetter' + str(j) + '_' + str(run)
        point = Point(data[j,0],data[j,1],
                    10)
        rate = data[j,2]
        # Convert rate to a value
        aval = np.log10(rate) + config['bvalue']*config["mmin"]

        mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, aval, config['bvalue'])
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
        source_list.append(point_source)

    mod_name = "Australia_Adaptive_K%i_b%.3f" % (smoother.config['k'], smoother.config['bvalue'])   
    nodes = list(map(obj_to_node, sorted(source_list)))
    source_model = Node("sourceModel", {"name": name}, nodes=nodes)
    with open(filename, 'wb') as f:
        nrml.write([source_model], f, '%s', xmlns = NAMESPACE)
                                      
# Set up paralell
proc = pypar.size()                # Number of processors as specified by mpirun                     
myid = pypar.rank()                # Id of of this process (myid in [0, proc-1])                     
node = pypar.get_processor_name()  # Host name on which current process is running                   
print 'I am proc %d of %d on node %s' % (myid, proc, node)
t0 = pypar.time()

config_params = params_from_shp(domains_shp, trt_ignore=['Interface', 'Active', 'Oceanic', 'Intraslab'])
for config in config_params:
    print config

#sys.exit()
# Read and clean the catalogue
parser = CsvCatalogueParser(ifile)
catalogue = parser.read_file(start_year=1900, end_year=2016)
# How many events in the catalogue?
print "The catalogue contains %g events" % catalogue.get_number_events()
neq = len(catalogue.data['magnitude'])
print "The catalogue contains %g events" % neq
# What is the geographical extent of the catalogue?
bbox = catalogue.get_bounding_box()
print "Catalogue ranges from %.4f E to %.4f E Longitude and %.4f N to %.4f N Latitude\n" % bbox
catalogue.sort_catalogue_chronologically()
index = np.logical_and(catalogue.data["magnitude"] > 1.5, catalogue.data["depth"] >= 0.0) 
#index = np.logical_and(catalogue.data["magnitude"] > 1.5, catalogue.data["magnitude"] < 4.0)
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
catalogue_depth_clean.get_number_events()

# Map configuration
map_config = {'min_lon': np.floor(100), 'max_lon': np.ceil(160),
              'min_lat': np.floor(-46), 'max_lat': np.ceil(-4), 'resolution':'c'}

grid_lims = [105., 160.0, 0.1, -47.0, -5.0, 0.1, 0., 20., 20.]

#try:
#    os.remove("Aus1_tmp.hdf5")
#except OSError:
#    pass
#config = {"bandwidth": 50,
#          "r_min": 1.0E-7, 
#          "bvalue": bvalue, "mmin": 3.0,
#          "learning_start": 1965, "learning_end": 2003,
#          "target_start": 2007, "target_end": 201
for i in range(0, len(config_params)*3, 1):
    if i % proc == myid:
        run = "%03d" % i
        print 'Run %s' % run
        completeness_table = config_params[i/3]['COMPLETENESS']
        if i % 3 == 0:
            bvalue = config_params[i/3]['BVAL_BEST']
        if i % 3 == 1:
            bvalue = config_params[i/3]['BVAL_UPPER']
        if i % 3 == 2:
            bvalue = config_params[i/3]['BVAL_LOWER']
        try:
            os.remove(("Aus1_tmp2%.3f_%s.hdf5" % (bvalue, run)))
        except OSError:
            pass
        mmin = completeness_table[0][1]
        print 'mmin', mmin
        config = {"k": 3,
                  "r_min": 1.0E-6, 
                  "bvalue": bvalue, "mmin": mmin,
                  "learning_start": 1900, "learning_end": 2016,
                  "target_start": 2018, "target_end": 2019} # using already optimised parameters
        ystart = completeness_table[-1][0]
        # Ensure we aren't training outside completeness model
        if ystart > config['learning_start']:
            config['learning_start'] = ystart

        run_smoothing(grid_lims, config, catalogue_depth_clean, completeness_table, map_config, run)

pypar.barrier()

if myid == 0:
    ss = int(pypar.time() - t0)
    h = ss / 3600
    m = (ss % 3600) / 60
    s = (ss % 3600) % 60
    print "--------------------------------------------------------"
    print 'P0: Total time (%i seconds): %s:%s:%s (hh:mm:ss)' % (ss,
                                                                string.zfill(h, 2),
                                                                string.zfill(m, 2),
                                                                string.zfill(s,2))
    print "--------------------------------------------------------"
pypar.finalize()
