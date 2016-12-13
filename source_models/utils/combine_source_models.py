# coding: utf-8
""" Combine multiple openquake input source models from different
files into one nrml format source model file
Jonathan Griffin, Geoscience Australia
April 2016
"""
# Python dependences
import os, time, glob
import numpy as np   # Numpy - Python's numerical library

# Input and Output Tools
from hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser  # Imports a source model from XML

# Load in the source models
source_folder = '../oq_inputs/source_group/'
source_model_list = []
for source_model_file in glob.glob(os.path.join(source_folder, '*.xml')):
#    print source_model_file
    parser = nrmlSourceModelParser(source_model_file)
    source_model_name = source_model_file.split('/')[-1]
#    print source_model_name
    source_model = parser.read_file(source_model_name) # You need to supply a name for the source model
    source_model_list.append(source_model)

# Add all sources to first source
source_list = []
source_model_index = 0
for source_model in source_model_list:
    source_index = 0
    for source in source_model.sources:
        # Take care of multiple source models with the same name here
        name = source_model.name + '_' + source.name + '_' + str(source_model_index) + \
               '_' + str(source_index)
        source.name = name
        print 'Adding source model ', source.name
        source_list.append(source)
        source_index += 1
    source_model_index += 1
#print source_list
source_model_list[0].sources = source_list
#print source_model_list[0].sources
# Write to nrml file 
source_model_list[0].serialise_to_nrml("Combined_Source_Model.xml")
