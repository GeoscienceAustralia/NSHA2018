
# coding: utf-8

# In[24]:

#get_ipython().magic(u'matplotlib inline')

# Python dependences
import os
import sys
import numpy as np   # Numpy - Python's numerical library
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # Matplotlib - Python's plotting library
from copy import deepcopy   # Python module for copying objects
import ogr
import shapefile
from shapely.geometry import Point, Polygon

# For running in parallel
import time
from time import localtime, strftime, gmtime
import string
import pypar

# Other smoothed seismicity functions
from utilities import params_from_shp

# Input and Output Tools
# Catalogue and sources 
from hmtk.parsers.catalogue import CsvCatalogueParser   # Reads an earthquake catalogue from CSV
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueWriter  # Writes an earthquake catalogue to CSV
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser  # Imports a source model from XML
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.sourcewriter import obj_to_node
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

# To build source model
from hmtk.sources.source_model import mtkSourceModel
from hmtk.sources.point_source import mtkPointSource
from openquake.hazardlib.scalerel.leonard2014 import Leonard2014_SCR
from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.tom import PoissonTOM
from openquake.hazardlib.sourcewriter import obj_to_node
from openquake.baselib.node import Node
from openquake.hazardlib import nrml
from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.mfd import TruncatedGRMFD
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
print "Everything Imported OK!"

domains_shp = '../zones/2018_mw/Domains_single_mc/shapefiles/Domains_NSHA18_MFD.shp'
#Importing catalogue
catalogue_filename = "../../catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv"
# Flag for whether to overwrite exiting .xml source model 
# files with the same b value and completeness combination.
# Shoudld normally set to True unless you are being really careful.
overwrite = True
#####################################
# Shouldn't need input below here
#####################################


def run_smoothing(grid_lims, smoothing_config, catalogue, completeness_table, map_config, 
                  run, overwrite=True):
    """Run all the smoothing
    """
    ystart = completeness_table[-1][0]
    yend = catalogue.end_year
    catalogue_comp = deepcopy(catalogue)
    # Ensuring that catalogue is cleaned of earthquakes outside of 
    # completeness period
    index = catalogue_comp.data['year']>=ystart
    catalogue_comp.purge_catalogue(index)

    completeness_string = 'comp'
    for ym in completeness_table:
        completeness_string += '_%i_%.1f' % (ym[0], ym[1])
    smoother_filename = 'Australia_Fixed_%i_%i_b%.3f_mmin_%.1f_0.1%s.csv' % (smoothing_config["BandWidth"], smoothing_config["Length_Limit"],
                                                                    bvalue, completeness_table[0][1], completeness_string)
    filename = smoother_filename[:-4] + '.xml'
    if os.path.exists(filename) and not overwrite:
        print '%s already created, not overwriting!' % filename
        return
    smoother = SmoothedSeismicity([105.,160.,0.1,-47.,-5,0.1,0.,20., 20.], bvalue = smoothing_config['bvalue'])
    print 'Running smoothing'
    smoothed_grid = smoother.run_analysis(catalogue_comp, smoothing_config, completeness_table=completeness_table)

    smoother.write_to_csv(smoother_filename)


    from openquake.hazardlib.nrml import SourceModelParser, write, NAMESPACE
    from openquake.baselib.node import Node
    from openquake.hazardlib import nrml
    from openquake.hazardlib.sourcewriter import obj_to_node
    # Build nrml input file of point sources
    source_list = []
    #i=0
    min_mag = 4.5
    max_mag = 7.8
    bval = bvalue # just define as 1 for time being
    # Read in data again to solve number fomatting issue in smoother.data
    # For some reason it just returns 0 for all a values
    try:
        data = np.genfromtxt(smoother_filename, delimiter = ',', skip_header = 1)
    except ValueError:
        print 'Something wrong with file %s' % smoother_filename
        sys.exit()
    tom = PoissonTOM(50) # Dummy temporal occurence model for building pt sources
    msr = Leonard2014_SCR()
    for j in range(len(data[:,4])):
    #    print smoother.data[j,:]
        identifier = 'FSS' + str(j) + '_' + str(run)
        name = 'Frankel' + str(j) + '_' + str(run)
        point = Point(data[j,0],data[j,1],
                    data[j,2])
        annual_rate = data[j,4]/(yend - ystart + 1)
        aval =  np.log10(annual_rate) + smoothing_config['bvalue']*completeness_table[0][1]
        mfd = TruncatedGRMFD(min_mag, max_mag, 0.1, aval, bval)
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

    nodes = list(map(obj_to_node, sorted(source_list)))
    source_model = Node("sourceModel", {"name": name}, nodes=nodes)
    with open(filename, 'wb') as f:
        nrml.write([source_model], f, '%s', xmlns = NAMESPACE)

    # Creating a basemap - input a cconfiguration and (if desired) a title
    title = 'Smoothed seismicity rate for learning \nperiod %i 2017, Mmin = %.1f' %(
         completeness_table[0][0], completeness_table[0][1])
    basemap1 = HMTKBaseMap(map_config, 'Smoothed seismicity rate')
    # Adding the smoothed grip to the basemap
    sym = (2., 3.,'cx')
    x,y = basemap1.m(smoother.data[:,0], smoother.data[:,1])
    basemap1.m.scatter(x, y, marker = 's', c = np.log10(smoother.data[:,4]), cmap = plt.cm.coolwarm, zorder=10, lw=0,
                       vmin=-6.5, vmax = 1.5 )
    basemap1.m.drawcoastlines(linewidth=1, zorder=50) # Add coastline on top
    basemap1.m.drawmeridians(np.arange(map_config['min_lat'], map_config['max_lat'], 5))
    basemap1.m.drawparallels(np.arange(map_config['min_lon'], map_config['max_lon'], 5))
    plt.colorbar(label='log10(Smoothed rate per cell)')
    plt.legend()
    figname = smoother_filename[:-4] + '_smoothed_rates_map.png'
    plt.savefig(figname)

# Set up paralell
proc = pypar.size()                # Number of processors as specified by mpirun                     
myid = pypar.rank()                # Id of of this process (myid in [0, proc-1])                     
node = pypar.get_processor_name()  # Host name on which current process is running                   
print 'I am proc %d of %d on node %s' % (myid, proc, node)
t0 = pypar.time()

parser = CsvCatalogueParser(catalogue_filename) # From .csv to hmtk

# Read and process the catalogue content in a variable called "catalogue"
catalogue = parser.read_file(start_year=1965, end_year=2016)

# How many events in the catalogue?
print "The catalogue contains %g events" % catalogue.get_number_events()

# What is the geographical extent of the catalogue?
bbox = catalogue.get_bounding_box()
print "Catalogue ranges from %.4f E to %.4f E Longitude and %.4f N to %.4f N Latitude\n" % bbox

catalogue.sort_catalogue_chronologically()
catalogue.data['magnitude']
index = catalogue.data['magnitude']>1.5

# Copying the catalogue and saving it under a new name "catalogue_clean"
catalogue_clean = deepcopy(catalogue)

# remove nan magnitudes
catalogue_clean.purge_catalogue(index)
catalogue_clean.sort_catalogue_chronologically()
catalogue_clean.data['magnitude']
catalogue_clean.data['year']
catalogue_clean.get_decimal_time()
catalogue_clean.data['longitude']

catalogue_depth_clean = deepcopy(catalogue_clean)
index = catalogue_depth_clean.data['depth']>=0.
catalogue_depth_clean.purge_catalogue(index)

grid_lims = [105., 160.0, 0.1, -47.0, -5.0, 0.1, 0., 20., 20.]

# Map configuration
map_config = {'min_lon': np.floor(105), 'max_lon': np.ceil(155),
              'min_lat': np.floor(-45), 'max_lat': np.ceil(-9), 'resolution':'c'}

magnitude_bin_width = 0.1  # In magnitude units
time_bin_width = 1.0 # In years

config_params = params_from_shp(domains_shp, trt_ignore=['Interface', 'Active', 'Oceanic', 'Intraslab'])
# Get unique combinations of config parameters to avoid
# running the same models multiple times (and running into file access issues)
config_combinations = []
for i in range(0, len(config_params), 1):
    completeness_table = config_params[i]['COMPLETENESS']
#    print completeness_table
    config_combo = [completeness_table, config_params[i]['BVAL_BEST']]
    print config_combo
    if not any(np.array_equal(l[0], config_combo[0]) and l[1]==config_combo[1] for l in config_combinations):
        config_combinations.append(config_combo)
    config_combo = [completeness_table, config_params[i]['BVAL_UPPER']]
    if not any(np.array_equal(l[0], config_combo[0]) and l[1]==config_combo[1] for l in config_combinations):
        config_combinations.append(config_combo)
    config_combo = [completeness_table, config_params[i]['BVAL_LOWER']]
    if not any(np.array_equal(l[0], config_combo[0]) and l[1]==config_combo[1] for l in config_combinations):
        config_combinations.append(config_combo)
print 'config_combinations'
print config_combinations
print len(config_combinations)

for i in range(0, len(config_combinations), 1):
    if i % proc == myid:
        run = "%03d" % i
        print 'Run %s' % run
        completeness_table = np.array(config_combinations[i][0])
        print 'Only using full period of catalogue'
        completeness_table = np.array([completeness_table[0]])
        bvalue = config_combinations[i][1]
#        if i % 3 == 0:
#            bvalue = config_params[i/3]['BVAL_BEST']
#        if i % 3 == 1:
#            bvalue = config_params[i/3]['BVAL_UPPER']
#        if i % 3 == 2:
#            bvalue = config_params[i/3]['BVAL_LOWER']
        mmin = completeness_table[0][1]
        print 'mmin', mmin
        config = {"BandWidth": 50.,
                  "Length_Limit": 3.,
                  "increment": False,
                  "bvalue":bvalue}
        ystart = completeness_table[-1][0]
        # Call the smoothing module
        run_smoothing(grid_lims, config, catalogue_depth_clean, completeness_table, map_config, run, overwrite)
            

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






