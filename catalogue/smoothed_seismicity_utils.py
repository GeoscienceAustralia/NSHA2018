#!usr/bin/env/python
"""
Smoothed seismicity and related utility functions
"""
import os
import h5py
import shapefile
import numpy as np
import collections
import matplotlib.pyplot as plt
from datetime import datetime
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser

class ProgressCounter(object):
    """
    Print out progress every N percent
    """
    def __init__(self, ntot, nval=10, timer=False):
        """
        Setup counter
        """
        self.vals = np.hstack([np.arange(0, 100, nval) * (ntot / 100),
                              ntot - 1])
        self.idx = np.cumsum(nval * np.ones_like(self.vals)) - nval
        self.cntr = 0
        self.target = self.vals[0]
        self.timer = timer

    def update(self, i):
        """
        Update - i.e. print if the counter hits one of the target cells
        """
        if i == self.target:
            if self.timer:
                print "%s %% - %s" % (self.idx[self.cntr],
                                      str(datetime.now()))
            else:
                print "%s %%" % self.idx[self.cntr]
            self.cntr += 1
            if self.cntr < len(self.vals):
                self.target = self.vals[self.cntr]

    def reset(self):
        """
        To re-use the counter again, re-set the counter back to 0
        """
        self.cntr = 0
        self.target = self.vals[0]


def get_catalogue_bounding_polygon(catalogue):
    '''
    Returns a polygon containing the bounding box of the catalogue
    '''
    upper_lon = np.max(catalogue.data['longitude'])
    upper_lat = np.max(catalogue.data['latitude'])
    lower_lon = np.min(catalogue.data['longitude'])
    lower_lat = np.min(catalogue.data['latitude'])

    return Polygon([Point(lower_lon, upper_lat), Point(upper_lon, upper_lat),
                    Point(upper_lon, lower_lat), Point(lower_lon, lower_lat)])


class Grid(collections.OrderedDict):
    @classmethod
    def make_from_list(cls, grid_limits):
        new = cls()
        new.update({'xmin': grid_limits[0],
                    'xmax': grid_limits[1],
                    'xspc': grid_limits[2],
                    'ymin': grid_limits[3],
                    'ymax': grid_limits[4],
                    'yspc': grid_limits[5],
                    'zmin': grid_limits[6],
                    'zmax': grid_limits[7],
                    'zspc': grid_limits[8]})
        return new

    @classmethod
    def make_from_catalogue(cls, catalogue, spacing, dilate):
        '''
        Defines the grid on the basis of the catalogue
        '''
        new = cls()
        cat_bbox = get_catalogue_bounding_polygon(catalogue)

        if dilate > 0:
            cat_bbox = cat_bbox.dilate(dilate)

        # Define Grid spacing
        new.update({'xmin': np.min(cat_bbox.lons),
                    'xmax': np.max(cat_bbox.lons),
                    'xspc': spacing,
                    'ymin': np.min(cat_bbox.lats),
                    'ymax': np.max(cat_bbox.lats),
                    'yspc': spacing,
                    'zmin': 0.,
                    'zmax': np.max(catalogue.data['depth']),
                    'zspc': np.max(catalogue.data['depth'])})

        if new['zmin'] == new['zmax'] == new['zspc'] == 0:
            new['zmax'] = new['zspc'] = 1

        return new

    def as_list(self):
        return [self['xmin'], self['xmax'], self['xspc'],
                self['ymin'], self['ymax'], self['yspc'],
                self['zmin'], self['zmax'], self['zspc']]

    def as_polygon(self):
        return Polygon([
            Point(self['xmin'], self['ymax']),
            Point(self['xmax'], self['ymax']),
            Point(self['xmax'], self['ymin']),
            Point(self['xmin'], self['ymin'])])

    def dilate(self, width):
        polygon = self.as_polygon().dilate(width)

        self.update({'xmin': np.min(polygon.lons),
                     'xmax': np.max(polygon.lons),
                     'ymin': np.min(polygon.lats),
                     'ymax': np.max(polygon.lats)})
        return self


def load_in_area_sources(filename):
    """
    Reads in an area source model in shapefile format and returns information
    as a set of ordered dictionaries
    """
    sf = shapefile.Reader(filename)
    fields = [fld[0] for fld in sf.fields[1:]]
    sources = []
    for rec in sf.shapeRecords():
        attribs = []
        for i, fld in enumerate(fields):
            attrib = rec.record[i].strip()
            if attrib:
                attribs.append((fields[i], attrib))
        attribs.append(("geometry", np.array(rec.shape.points)))
        sources.append(collections.OrderedDict(attribs))
    return sources

def build_spatio_temporal_completeness_vectors(catalogue, sources,
        completeness_file, mmin, default_bval=1.0):
    """
    
    """
    # Associated completeness tables to sources
    comp_fle = h5py.File(completeness_file)
    for src in sources:
        # Add the completeness table to the sources
        if src["id"] in comp_file.keys():
            src_comp_table = comp_file[src["id"]][:]
            src["completeness"] = src_comp_table
        else:
            print "No completeness table for source %s" % src["id"]
    # Now process the catalogue to associate each event with a source
    event_sources = [None for i in range(catalogue.get_number_events())]
    weights = np.zeros(catalogue.get_number_events())
    lonlats = np.column_stack([catalogue.data["longitude"],
                               catalogue.data["latitude"]])
    for src in sources:
        # Find earthquakes in source
        src_path = Path(src["geometry"])
        inpoly = np.where(src_path.contains(lonlats))[0]
        if len(inpoly):
            for idx in inpoly:
                event_sources[idx] = src["id"]
            if "completeness" in src.keys():
                # Get completeness vectors
                cyear = np.hstack([np.inf, src["completeness"][:, 0]])
                cmag = src["completeness"][:, 1]
                for i in range(len(cyear) - 1):
                    idx = np.logical_and(
                        catalogue.data["year"] < cyear[i],
                        catalogue.data["year"] >= cyear[i + 1]
                        )
                    idx = np.logical_and(
                        idx,
                        catalogue.data["magnitude"] >= cmag[i]
                        )
                    if "b_val" in src.keys():
                        bval = src["b_val"]
                    else:
                        bval = default_bval
                    weights[idx] = 10.0 ** (bval * (cmag[i] - mmin))
    return weights

                



        
            


