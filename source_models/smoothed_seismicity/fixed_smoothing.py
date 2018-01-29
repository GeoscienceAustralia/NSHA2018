
# coding: utf-8

# In[24]:

#get_ipython().magic(u'matplotlib inline')

# Python dependences
import os
import sys
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
#from nrml.models import HypocentralDepth
#nrml.models import HypocentralDepth
print "Everything Imported OK!"


# In[25]:
bvalue = float(sys.argv[1])
print 'b value', bvalue

#Importing catalogue
catalogue_filename = "../../catalogue/data/AUSTCAT_V0.12_hmtk_declustered.csv"
#catalogue_filename = "../../catalogue/data/AUSTCAT_V0.12_hmtk_mx_orig.csv"
parser = CsvCatalogueParser(catalogue_filename) # From .csv to hmtk

# Read and process the catalogue content in a variable called "catalogue"
catalogue = parser.read_file(start_year=1965, end_year=2010)

# How many events in the catalogue?
print "The catalogue contains %g events" % catalogue.get_number_events()

# What is the geographical extent of the catalogue?
bbox = catalogue.get_bounding_box()
print "Catalogue ranges from %.4f E to %.4f E Longitude and %.4f N to %.4f N Latitude\n" % bbox


# In[26]:

catalogue.sort_catalogue_chronologically()
catalogue.data['magnitude']
index = catalogue.data['magnitude']>1.5
index
#min(catalogue.data['depth'])


# In[27]:

# Copying the catalogue and saving it under a new name "catalogue_clean"
catalogue_clean = deepcopy(catalogue)

# remove nan magnitudes
catalogue_clean.purge_catalogue(index)
catalogue_clean.sort_catalogue_chronologically()
catalogue_clean.data['magnitude']
catalogue_clean.data['year']
catalogue_clean.get_decimal_time()
catalogue_clean.data['longitude']


# In[28]:

catalogue_depth_clean = deepcopy(catalogue_clean)
index = catalogue_depth_clean.data['depth']>=0.
catalogue_depth_clean.purge_catalogue(index)


# In[29]:

catalogue_clean.get_number_events()


# In[30]:

#source_model_file = "../zones/2012_mw_ge_4.0/NSHA13_Background/input/best/NSHA13_BACKGROUND_best.xml"
source_model_file = 'Aus_cont_testzone.xml'
parser = nrmlSourceModelParser(source_model_file)

# Parse the seismic sources and save them into a variable called "source_model"
source_model = parser.read_file("Aus Source Model 1") # You need to supply a name for the source model


# In[31]:

#plot_magnitude_time_scatter(catalogue_clean, plot_error=False)


# In[32]:

# Map configuration
llon, ulon, llat, ulat = catalogue_clean.get_bounding_box()
#map_config = {'min_lon': np.floor(llon), 'max_lon': np.ceil(ulon),
 #             'min_lat': np.floor(llat), 'max_lat': np.ceil(ulat), 'resolution':'c'}
map_config = {'min_lon': np.floor(105), 'max_lon': np.ceil(155),
              'min_lat': np.floor(-45), 'max_lat': np.ceil(-9), 'resolution':'c'}
# Creating a basemap - input a cconfiguration and (if desired) a title
#basemap1 = HMTKBaseMap(map_config, 'Earthquake Catalogue')

# Adding the seismic sources
#basemap1.add_source_model(source_model, area_border='r-', border_width=1.5, alpha=0.5)


# In[33]:

# Select catalogue from within sourcezone
selector1 = CatalogueSelector(catalogue_depth_clean, create_copy=True)
for source in source_model.sources:
    source.select_catalogue(selector1)
    
    llon, ulon, llat, ulat = source.catalogue.get_bounding_box()
    print llon, ulon, llat, ulat
    # Map the Source
    #src_basemap = HMTKBaseMap(map_config, "Source: {:s}".format(source.name))
    #print "Source ID: %s  Source Name: %s   Number of Events: %g" % (source.id, source.name,
    #                                                                 source.catalogue.get_number_events())
    # Add on the catalogue
    #src_basemap.add_catalogue(source.catalogue, overlay=False)


    


# In[34]:

# Map configuration
llon, ulon, llat, ulat = catalogue_clean.get_bounding_box()
#map_config = {'min_lon': np.floor(llon), 'max_lon': np.ceil(ulon),
#              'min_lat': np.floor(llat), 'max_lat': np.ceil(ulat), 'resolution':'c'}
map_config = {'min_lon': np.floor(100), 'max_lon': np.ceil(160),
              'min_lat': np.floor(-45), 'max_lat': np.ceil(-4), 'resolution':'c'}
# Creating a basemap - input a cconfiguration and (if desired) a title
#basemap1 = HMTKBaseMap(map_config, 'Earthquake Catalogue')


# Adding the catalogue to the basemap
# In this case we will 'close' the figure after rendering, we do this by setting 'overlay=False'
# This is also the default option
# If we wanted to add another layer on top, we would set the overlay to True
#basemap1.add_catalogue(catalogue_depth_clean, overlay=False)


# In[35]:

magnitude_bin_width = 0.1  # In magnitude units
time_bin_width = 1.0 # In years
#plot_magnitude_time_density(source.catalogue, magnitude_bin_width, time_bin_width)


# In[36]:

# Shows depth histogram every 5 km  
#plot_depth_histogram(source.catalogue, 5., normalisation=True)


# In[37]:

completeness_table_a = np.array([[1965., 4.0]])
#completeness_table_a = np.array([[1978., 3.5],
#                                 [1964., 4.0],
#                                 [1900., 5.5]])
#completeness_table_a = np.array([[1900., 5.0]])
#plot_magnitude_time_density(source.catalogue, 0.1, 1.0,
#                            completeness=completeness_table_a)


# In[38]:

#from hmtk.seismicity.occurrence.weichert import Weichert

#recurrence_estimator = Weichert()

#recurrence_config = {"magnitude_interval": 0.1}

#bval, sigma_b, aval, sigma_a = recurrence_estimator.calculate(source.catalogue,
#                                                              recurrence_config,
#                                                              completeness_table_a)

#print "a = %.3f (+/- %.3f),  b = %.3f (+/-%.3f)" % (aval, sigma_a, bval, sigma_b)


# In[39]:

#from hmtk.plotting.seismicity.occurrence.recurrence_plot import plot_recurrence_model
#from openquake.hazardlib.mfd import TruncatedGRMFD
#mfd0 = TruncatedGRMFD(4.5, 7.2, 0.1, aval, bval)#
#plot_recurrence_model(mfd0, source.catalogue, completeness_table_a, 0.1)


# In[ ]:

from hmtk.seismicity.smoothing.smoothed_seismicity import SmoothedSeismicity
smoothing_config = {"BandWidth": 50.,
                    "Length_Limit": 3.,
                    "increment": 0.1}
#bvalue = 0.819
#bvalue = 0.835
#upper
#bvalue = 0.747
#bvalue= 0.727
#bvalue= 1.062
#lower
#bvalue = 0.892
#bvalue = 0.944
#bvalue = 1.355
smoother = SmoothedSeismicity([100.,160.,0.1,-45.,-5,0.1,0.,20., 20.], bvalue = bvalue)
#smoothed_grid = smoother.run_analysis(source_model.sources[0].catalogue, smoothing_config, completeness_table=completeness_table_a)
print 'Running smoothing'
smoothed_grid = smoother.run_analysis(source.catalogue, smoothing_config, completeness_table=completeness_table_a)
smoother_filename = 'smoothed_%i_%i_mmin_%.1f_%.3f_0.1.csv' % (smoothing_config["BandWidth"], smoothing_config["Length_Limit"],
                                                               completeness_table_a[0][-1], bvalue)
smoother.write_to_csv(smoother_filename)


# In[ ]:

#smoother_filename = 'smoothed_%i_%i_mmin_%.1f_0.1.csv' % \
##                        (smoothing_config["BandWidth"], smoothing_config["Length_Limit"],
  #                      completeness_table_a[0][-1])
#smoother.write_to_csv(smoother_filename)


# In[ ]:

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
data = np.genfromtxt(smoother_filename, delimiter = ',', skip_header = 1)
#print max(data[:,4])
#print data[:,4]
#print len(data[:,4])
tom = PoissonTOM(50) # Dummy temporal occurence model for building pt sources
msr = Leonard2014_SCR()
for j in range(len(data[:,4])):
#    print smoother.data[j,:]
    identifier = 'FSS' + str(j)
    name = 'Frankel' + str(j)
    point = Point(data[j,0],data[j,1],
                data[j,2])
    rate = data[j,4]
    aval = np.log10(rate)
   # aval = rate # trying this based on some testing
#    aval = np.log10(rate) #+ bval*completeness_table_a[0][1]
   # print aval
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
#    i+=1
#    if j==1000:
#        break

filename = "smoothed_frankel_50_3_mmin_%.1f_b%.3f_0.1.xml" % (completeness_table_a[0][-1], bvalue)
mod_name = 'smoothed_frankel_50_3_mmin_%.1f_b%.3f_0.1' % (completeness_table_a[0][-1], bvalue)    
nodes = list(map(obj_to_node, sorted(source_list)))
source_model = Node("sourceModel", {"name": name}, nodes=nodes)
with open(filename, 'wb') as f:
    nrml.write([source_model], f, '%s', xmlns = NAMESPACE)

#source_model = mtkSourceModel(identifier=0, name='Frankel_50_3',
#                              sources = source_list)
#source_model.serialise_to_nrml('smoothed_frankel_50_3_mmin_%.1f_b%.3f_0.1_full.xml' % (completeness_table_a[0][-1], bvalue))


# In[ ]:

# Map configuration
#llon, ulon, llat, ulat = source.catalogue.get_bounding_box()
#map_config = {'min_lon': np.floor(llon), 'max_lon': np.ceil(ulon),
#              'min_lat': np.floor(llat), 'max_lat': np.ceil(ulat), 'resolution':'c'}
#map_config = {'min_lon': np.floor(105), 'max_lon': np.ceil(155),
#              'min_lat': np.floor(-45), 'max_lat': np.ceil(-9), 'resolution':'c'}
# Creating a basemap - input a cconfiguration and (if desired) a title
#basemap1 = HMTKBaseMap(map_config, 'Smoothed seismicity rate')
#basemap1.m.drawmeridians(np.arange(llat, ulat, 5))
#basemap1.m.drawparallels(np.arange(llon, ulon, 5))
#print smoother.data[:,0]
#print smoother.data[:,1]
# Adding the smoothed grip to the basemap
#sym = (2., 3.,'cx')
#x,y = basemap1.m(smoother.data[:,0], smoother.data[:,1])
#print data[:,4]
#basemap1.m.scatter(x, y, marker = 's', c = np.log10(smoother.data[:,4]), cmap = plt.cm.coolwarm, zorder=10, lw=0,
#                   vmin=-6.5, vmax = 1.5 )
#basemap1.m.scatter(x, y, marker = 's', c = np.arange(-7.5, -0.5, 0.1), cmap = plt.cm.coolwarm, zorder=10, lw=0)
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
#plt.savefig('smoothed_%i_%i_%.1f_b%.3f_0.1.png' %                         (smoothing_config["BandWidth"], smoothing_config["Length_Limit"],                         completeness_table_a[0][-1], bvalue))


# In[ ]:



