#!/usr/bin/env/python
"""
Trial implementation of Helmstetter & Werner (2012)
From Graeme Weatherill, GEM Foundation, with minor changes
by Jonathan Griffin, Geoscience Australia
"""
import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
from  matplotlib import pyplot as plt
import collections
from matplotlib.path import Path
from copy import deepcopy
from datetime import datetime 
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo import Point, Line, Polygon
from scipy.special import erf
from scipy.optimize import minimize
from scipy.stats import norm, truncnorm
from scipy.misc import factorial
from math import fabs, floor, log10, sqrt, atan2, pi, degrees, radians
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser
from smoothed_seismicity_utils import (ProgressCounter, Grid, 
                                       get_catalogue_bounding_polygon)

SAM_COMP_TABLE = np.array([[1995.0, 3.5],
                           [1985.0, 4.5],
                           [1964.0, 5.0],
                           [1950.0, 6.0],
                           [1920.0, 6.5],
                           [1910.0, 7.0],
                           [1900.0, 7.5]])

def _cnn(p, a):
    return p[0] + a * p[1]

def _time_constraint(p, a, k, data):
    """
    Defines a time-based constraint on a set
    """
    k = k if (k >= 1) else 1.0
    time_set = data["time"][data["nearest"][1:(k + 1)]]
    if np.all(time_set < 0.0):
        return p[0] - np.min(time_set)
    else:
        return p[0] - np.max(np.fabs(time_set))

    #return p[0] - np.max(np.fabs(data["time"][data["nearest"][1:(k + 1)]]))
    #return np.min(data["time"][data["nearest"][1:(k + 1)]) + p[0]

def _space_constraint(p, a, k, data):
    """
    Defines a spatial constraint on the set

    """
    k = k if (k >= 1) else 1.0
    # Get nearest k events
    return p[1] - np.max(data["distance"][data["nearest"][1:(k + 1)]])


def optimise_cnn(initial_p, a, k, data):
    """
    Optimise the parameters of the coupled nearest neighbour algorithm
    """
    consts = [
        {"type": 'ineq', 'fun': _time_constraint, 'args': [a, k, data]},
        {"type": 'ineq', 'fun': _space_constraint, 'args': [a, k, data]}]
    bnds = [(0., None), (0.0, None)]
    result = minimize(_cnn, initial_p, args=[a],
                      method="SLSQP",
                      bounds=bnds,
                      constraints=consts)
    return result.x
                        

def _get_length_scales(eqlon, eqlat, bandwidth, rad_earth=6371.0):
    '''
    Returns the bandwidth for the event in terms of decimal degrees
    :param float eqlon:
        Earthquake longitude
    :param float eqlat:
        Earthquake latitude
    :param float bandwidth:
        Bandwidth of smoothing kernel (km)
    '''
    xlen_scale = degrees(bandwidth / (rad_earth * np.cos(radians(eqlat))))
    ylen_scale = (bandwidth / (2. * pi * rad_earth)) * 360.
    return xlen_scale * sqrt(2.), ylen_scale * sqrt(2.)


def _get_distances(eqlon, eqlat, clons, clats, spcx, spcy):
    '''
    Returns the distances (xeq - xcell) and (yeq - ycell) taking into account
    the crossing of the International Dateline
    '''
    if not spcy:
        spcy = spcx
    lx = clons - (spcx / 2.)
    ux = clons + (spcx / 2.)
    # Get longitude length scales
    dx1 = eqlon - lx
    dx2 = eqlon - ux
    # Address cases of lower bound across dateline and earthquake is in
    # eastern hemisphere
    idx = np.logical_and(np.fabs(dx1) > 180., dx1 > 0.)
    #dx1[idx] = -(360. - dx1[idx]) 
    dx1[idx] = eqlon - (360 + lx[idx])
    dx2[idx] = eqlon - (360 + ux[idx])
    
    #dx2[idx] = dx1[idx] + spcx
    # Address cases of lower abound across dateline and earthquake is in 
    # western hemisphere
    idx = np.logical_and(np.fabs(dx1) > 180., dx1 < 0.)
    dx1[idx] = eqlon - (-360. + lx[idx]) 
    dx2[idx] = eqlon - (-360. + ux[idx]) 
    #dx2[idx] = dx1[idx] - spcx
    
    # For latitudes this is not a problem - just take the values
    dy1 = eqlat - (clats - (spcy / 2.))
    dy2 = eqlat - (clats + (spcy / 2.))

    # For depth return differences in kilometres 
    #dz1 = eqdep - (cdeps - (spcz / 2.))
    #dz2 = eqdep - (cdeps + (spcz / 2.))

    return dx1, dx2, dy1, dy2


class HelmstetterWerner2012(object):
    """
    This method implements the smoothed seismicity method of Helmstetter &
    Werner (2012)

    :param list grid_limits:
        Defines the limits of the smoothing grid: [llon, ulon, dlon, llat,
        ulat, dlat, udepth, ldepth, ddepth]
    :param dict config:
        Dictionary of configuration parameters
        Configuration parameters: 
        - h0 - Time bandwidth (days)
          d0 - Distance bandwidth (km)
          rmin - Baseline rate (rate / cell / day)
          bvalue - b-value
          learning_start - Start of learning catalogue (year)
          learning_end - End of learning catalogue (year)
          target_start - Start of target catalogue
          target_end - End of target catalogue
          grid_limits - Bounding box of grid limits
          catalogue - Earthquake catalogue
    :param numpy.ndarry lon_bins:
        Longitude bin edges
    :param numpy.ndarry lat_bins:
        Latitude bin edges
    :param numpy.ndarry depth_bins:
        Depth bin edges
    :param fle:
        h5py.File object holding distance info
    :param learning_weights:
        completeness weights of the learning catalogue
    :param target_weights:
        completeness weights of the target catalogue
    :param catalogue:
        Earthquake catalogue as instance of the class
        hmtk.seismicity.catalogue.Catalogue
    :param learning_catalogue:
        Learning portion of the earthquake catalogue as instance of the class
        hmtk.seismicity.catalogue.Catalogue
    :param target_catalogue:
        Target portion of the earthquake catalogue as instance of the class
        hmtk.seismicity.catalogue.Catalogue
    :param numpy.ndarray rates:
        smoothed seismicity rates
    """
    def __init__(self, grid_limits, config, catalogue,
            storage_file="tmp.hdf5"):
        """
        Instantiate
        :param str storage_file:
            Path to file for storing results
        """
        self.grid_limits = Grid.make_from_list(grid_limits)
        self.config = config
        self.lon_bins = np.arange(
            self.grid_limits["xmin"],
            self.grid_limits["xmax"] + self.grid_limits["xspc"],
            self.grid_limits["xspc"])
        self.lat_bins = np.arange(
            self.grid_limits["ymin"],
            self.grid_limits["ymax"] + self.grid_limits["yspc"],
            self.grid_limits["yspc"])
        self.depth_bins = np.arange(
            self.grid_limits["zmin"],
            self.grid_limits["zmax"] + self.grid_limits["zspc"],
            self.grid_limits["zspc"])
        self._build_grid()
        #if not os.path.exists(storage_file):
        #    os.path.mkdir(temp_dir)
        self.fle = h5py.File(storage_file)
        self.learning_weights = None
        self.target_weights = None
        self.catalogue = catalogue
        self.catalogue.data["dtime"] = catalogue.get_decimal_time()
        self.learning_catalogue = None
        self.target_catalogue = None
        self._split_catalogue()
        self.rates = None
        self.opt_counter = 0

    def _split_catalogue(self):
        """
        Given the configuration split the catalogue into a learning and
        a target catalogue
        """
        idx_l = np.logical_and(
            self.catalogue.data["year"] >= self.config["learning_start"],
            self.catalogue.data["year"] <= self.config["learning_end"])
        
        idx_l = np.logical_and(
            self.catalogue.data["magnitude"] >= self.config["mmin"],
            idx_l)

        idx_t = np.logical_and(
            self.catalogue.data["year"] >= self.config["target_start"],
            self.catalogue.data["year"] <= self.config["target_end"])
        idx_t = np.logical_and(
            self.catalogue.data["magnitude"] >= self.config["mmin"],
            idx_t)
        self.learning_catalogue = deepcopy(self.catalogue)
        self.learning_catalogue.purge_catalogue(idx_l)
        self.learning_catalogue.start_year = self.config["learning_start"]
        self.learning_catalogue.end_year = self.config["learning_end"]
        self.target_catalogue = deepcopy(self.catalogue)
        self.target_catalogue.purge_catalogue(idx_t)
        self.target_catalogue.start_year = self.config["target_start"]
        self.target_catalogue.end_year = self.config["target_end"]

    def _build_grid(self):
        """
        Function to build grid
        """
        g_x, g_y = np.meshgrid(
            self.lon_bins[:-1] + (self.grid_limits["xspc"] / 2.),
            self.lat_bins[:-1] + (self.grid_limits["yspc"] / 2.))
        n_x, n_y = g_x.shape
        ngp = n_x * n_y
        g_x = np.reshape(g_x, [ngp, 1]).flatten()
        g_y = np.reshape(g_y, [ngp, 1]).flatten()
        self.grid = np.column_stack([g_x, g_y])
        self.ngpts = self.grid.shape[0]

    def _get_catalogue_completeness_weights(self, completeness_table):
        """
        Returns the catalogue completeness weights based on the formulation
        w_i = 10 ^ (b * (m_c(t, r) - mmin))
        :param completeness_table:
            Completeness as a table of [[year, mag]]
        """
        cyear = np.hstack([np.inf, completeness_table[:, 0]])
        cmag = completeness_table[:, 1]
        self.learning_weights = np.zeros(
            self.learning_catalogue.get_number_events())
        self.target_weights = np.zeros(
            self.target_catalogue.get_number_events())
        # Define weights for complete magnitudes
        for i in range(len(cyear) - 1):
            idx = np.logical_and(
                self.learning_catalogue.data["year"] < cyear[i],
                self.learning_catalogue.data["year"] >= cyear[i + 1])
            idx = np.logical_and(
                idx,
                self.learning_catalogue.data["magnitude"] >= cmag[i])
            self.learning_weights[idx] = 10.0 ** (
                self.config["bvalue"] * (cmag[i] - self.config["mmin"])
                )
        # Purge learning catalogue of incomplete events
        self.learning_catalogue.purge_catalogue(self.learning_weights > 0.0)
        self.learning_weights = self.learning_weights[
            self.learning_weights > 0.0]
        # Same for target catalogue
        for i in range(len(cyear) - 1):
            idx = np.logical_and(
                self.target_catalogue.data["year"] < cyear[i],
                self.target_catalogue.data["year"] >= cyear[i + 1])
            idx = np.logical_and(
                idx,
                self.target_catalogue.data["magnitude"] >= cmag[i])
            self.target_weights[idx] = 10.0 ** (
                self.config["bvalue"] * (cmag[i] - self.config["mmin"])
                )
        self.target_catalogue.purge_catalogue(self.target_weights > 0.0)
        self.target_weights = self.target_weights[self.target_weights > 0.0]

    def _get_spatio_temporal_completeness_weights(self, completeness_file):
        """
        Returns the completeness weights assuming a spatio-temporal
        completeness model defined as an hdf5 file with each completeness
        region given as a group containing the datasets "Completeness" and
        "Outline" containing the completeness table and the outline of the
        zones respectively
        """
        cfle = h5py.File(completeness_file)
        self.learning_weights = np.zeros(
            self.learning_catalogue.get_number_events())
        self.target_weights = np.zeros(
            self.target_catalogue.get_number_events())
        lonlat = np.column_stack([
            self.learning_catalogue.data["longitude"],
            self.learning_catalogue.data["latitude"]])
        for key in cfle.keys():
            outline = Path(cfle["/".join([key, "Outline"])][:])
            completeness_table = cfle["/".join([key, "Completeness"])][:]
            cyear = np.hstack([np.inf, completeness_table[:, 0]])
            cmag = completeness_table[:, 1]
            # Define weights for complete magnitudes
            for i in range(len(cyear) - 1):
                idx = np.logical_and(
                    self.learning_catalogue.data["year"] < cyear[i],
                    self.learning_catalogue.data["year"] >= cyear[i + 1])
                idx = np.logical_and(
                    idx,
                    self.learning_catalogue.data["magnitude"] >= cmag[i])
                idx = np.logical_and(idx, outline.contains_points(lonlat))
                self.learning_weights[idx] = 10.0 ** (
                    self.config["bvalue"] * (cmag[i] - self.config["mmin"])
                    )
        # Purge learning catalogue of incomplete events
        self.learning_catalogue.purge_catalogue(self.learning_weights > 0.0)
        self.learning_weights = self.learning_weights[
            self.learning_weights > 0.0]
        # Same for target catalogue
        lonlat = np.column_stack([
            self.target_catalogue.data["longitude"],
            self.target_catalogue.data["latitude"]])
        for key in cfle.keys():
            outline = Path(cfle["/".join([key, "Outline"])][:])
            completeness_table = cfle["/".join([key, "Completeness"])][:]
            cyear = np.hstack([np.inf, completeness_table[:, 0]])
            cmag = completeness_table[:, 1]
            for i in range(len(cyear) - 1):
                idx = np.logical_and(
                    self.target_catalogue.data["year"] < cyear[i],
                    self.target_catalogue.data["year"] >= cyear[i + 1])
                idx = np.logical_and(
                    idx,
                    self.target_catalogue.data["magnitude"] >= cmag[i])
                idx = np.logical_and(idx, outline.contains_points(lonlat))
                self.target_weights[idx] = 10.0 ** (
                    self.config["bvalue"] * (cmag[i] - self.config["mmin"])
                    )
        self.target_catalogue.purge_catalogue(self.target_weights > 0.0)
        self.target_weights = self.target_weights[self.target_weights > 0.0]
        cfle.close()

    def build_time_distance_arrays(self):
        """
        Determines the arrays of distance between events and from
        events to target sites
        """
        neq = self.learning_catalogue.get_number_events()
        print "Building Time arrays"
        time_dset = self.fle.create_dataset("time", (neq, neq), dtype="f")
        #nearest_tset = self.fle.create_dataset("nearest_time", (neq, neq),
        #                                       dtype="i")
        counter1 = ProgressCounter(neq, timer=True)
        for i in range(neq - 1):
            counter1.update(i)
            # Get time difference (in DAYS)
            time_kernel = (365.25 * 
                (self.learning_catalogue.data["dtime"] -
                self.learning_catalogue.data["dtime"][i]))
            time_dset[i, :] = time_kernel
            #nearest_tset[i, :] = np.argsort(time_kernel)

        print "Building Distance Arrays"
        distance_dset = self.fle.create_dataset("distance",
                                                (neq, neq),
                                                dtype="f")
        nearest_dset = self.fle.create_dataset("d_nearest",
                                               (neq, neq),
                                               dtype="i")
        counter1.reset()
        for i in range(neq - 1):
            counter1.update(i)
            epi_distance = geodetic.geodetic_distance(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                self.learning_catalogue.data["longitude"],
                self.learning_catalogue.data["latitude"])
            distance_dset[i, :] = epi_distance
            nearest_dset[i, :] = np.argsort(epi_distance)

    def build_catalogue_2_grid_array(self, etime, ndays):
        """
        In the hdf5 file store, for each earthquake, a Ngd by 4 array of
        dx1, dx2, dy1, dy2
        """
        d_t = SECONDS_PER_DAY * float(self.config["ndays"]) / SECONDS_PER_YEAR
        neq = self.learning_catalogue.get_number_events()
        print "Building event to grid distance arrays"
        print str(datetime.now())
        gdist_dset = self.fle.create_dataset("grid_distance",
            (self.ngpts, 4, neq), chunks=(self.ngpts, 4, 1), dtype="f")
        print "Building event to grid time arrays"
        print str(datetime.now())
        target_times = np.arange(etime, etime + (11.0 * d_t), d_t)
        tdist_dset = self.fle.create_dataset("grid_time",
            (neq, len(target_times)),
            dtype="f")
        print "Processing ..."
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        for i in range(neq):
            counter1.update(i)
            #print i
            # Get distances
            dx1, dx2, dy1, dy2 = _get_distances(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                self.grid[:, 0], 
                self.grid[:, 1],
                self.grid_limits["xspc"],
                self.grid_limits["yspc"])
            #print "Calc distance %s" % str(datetime.now())
            gdist_dset[:, :, i] = np.column_stack([dx1, dx2, dy1, dy2])
            #print "Write distance %s" % str(datetime.now())
            # Get time difference (in decimal days)
            #print "Calc time %s" % str(datetime.now())
            dtime = target_times - self.learning_catalogue.data["dtime"][i]
            #print "Write time %s" % str(datetime.now())
            tdist_dset[i, :] = dtime
        counter1.reset()
        print "done!"

    def optimise_bandwidths(self, maxiter=50):
        """
        Using the epicentre to epicentre rates
        """
        neq = self.learning_catalogue.get_number_events()
        inv_n = 1.0 / float(neq)
        h_i = self.config["h0"] * np.ones(neq)
        d_i = self.config["d0"] * np.ones(neq)
        distance = -(self.fle["distance"][:] ** 2.)
        time_diff = -(self.fle["time"][:] ** 2.)

        diff_lr = 1.0
        iter_count = 0.0
        cur_lr = -np.inf
        not_converged = True

        while not_converged:
            print "Optimising bandwidth - iteration %g" % (iter_count + 1)
            rates = self.config["r_min"] + np.zeros(neq)

            counter1 = ProgressCounter(neq, 5.0, timer=True)

            for i in range(self.learning_catalogue.get_number_events()):
                counter1.update(i)
                rates[:i] += (
                    (2. / (h_i[i] * (d_i[i] ** 2.))) *
                    (np.exp(time_diff[:i, i]) ** (1. / (h_i[i] ** 2.0))) *
                    (np.exp(distance[:i, i]) ** (1. / (d_i[i] ** 2.0)))
                    )
            counter1.reset()
            new_lr = np.sum(np.log(rates))
            print "Iteration %g - Likelihood %.8f" % (iter_count, new_lr)
            if new_lr > cur_lr:
                cur_lr = np.copy(new_lr)
                # Get the geometric mean of the rates
                geo_mean_rate = np.sqrt(
                    np.exp(np.mean(np.log(rates))) / rates
                    )
                h_i *= geo_mean_rate
                d_i *= geo_mean_rate
                iter_count += 1
            else:
                not_converged = False
            if iter_count > maxiter:
                print "Failed to converge!"
                not_converged = False
        return h_i, d_i

    def run_smoothing(self, r_min, h_i, d_i):
        """
        Runs the smoothing algorithm
        """
        ntime = self.fle["grid_time"].shape[1]
        self.rates = r_min * np.ones([self.ngpts, ntime])
        time_diff = -((self.fle["grid_time"][:]).astype("float64")) ** 2.
        neq = self.learning_catalogue.get_number_events()
        counter1 = ProgressCounter(neq, 5.0, timer=True)

        for i in range(self.learning_catalogue.get_number_events()):
            counter1.update(i)
            weight = self.learning_weights[i]
            dgeo = (self.fle["grid_distance"][:, :, i]).astype("float64")
            # time kernel
            time_kernel = (2.0 / h_i[i]) * np.exp(
                time_diff[i, :] / (h_i[i] ** 2.)
                )
            #print self.learning_catalogue.data["dtime"][i],\
            #    self.learning_catalogue.data["longitude"][i],\
            #    self.learning_catalogue.data["latitude"][i],\
            #    self.learning_catalogue.data["magnitude"][i],\
            #    time_kernel
            # Distance kernel
            # Get length scales based on bandwidth
            xls, yls = _get_length_scales(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                d_i[i])
            dx1 = dgeo[:, 0]
            dx2 = dgeo[:, 1]
            dy1 = dgeo[:, 2]
            dy2 = dgeo[:, 3]
            valid_long = np.logical_and(dx1 <= (5.92 * xls),
                                        dx2 >= (-5.92 * xls))
            valid_lat = np.logical_and(dy1 <= (5.92 * yls),
                                       dy2 >= (-5.92 * yls))
            idx = np.logical_and(valid_long, valid_lat)
            #print np.sum(valid_long), np.sum(valid_lat), np.sum(idx)
            distance_kernel = 0.25  * (
                erf(dx2[idx] / xls) - erf(dx1[idx] / xls)) *\
                (erf(dy2[idx] / yls) - erf(dy1[idx] / yls))
            for j in range(ntime):
                self.rates[idx, j] += (weight *
                                       time_kernel[j] * distance_kernel)
            #if np.any(distance_kernel > 1.0):
            #    print i, weight, time_kernel, np.max(distance_kernel)
            #if np.any(time_kernel > 1.0):
            #    print " ".join([str(i), str(h_i[i]), str(d_i[i]), str(weight),
            #                    ",".join([str(val) for val in time_kernel])])
        self.rates = np.median(self.rates, axis=1)
        # Convert to annual rates
        self.rates *= (365.25 / float(self.config["ndays"]))


class HelmstetterNearestNeighbour(HelmstetterWerner2012):
    """
    Adaptation of the Helmstetter & Werner smoothing algorithm based on
    a nearest neighbour approach
    """
    def optimise_bandwidths(self):
        """
        Nearest neighbour optimisation
        """
        initial_p = [10., 10.]
        # Build data from file
        neq = self.learning_catalogue.get_number_events()
        h_i = np.zeros(neq)
        d_i = np.zeros(neq)
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        for i in range(neq):
            counter1.update(i)
            data = {"time": self.fle["time"][i, :],
                    "distance": self.fle["distance"][i, :],
                    "nearest": self.fle["d_nearest"][i, :]}
            h_i[i], d_i[i] = optimise_cnn(initial_p,
                                          self.config["a"],
                                          self.config["k"],
                                          data)
        d_i[d_i < 0.5] = 0.5
        h_i[h_i < 10.] = np.median(h_i)
        return h_i, d_i

    def optimum_smoothing(self, params):
        """
        Runs the Helmstetter Adaptive Kernel Smoothing method optimized 
        """
        self.config["k"] = int(params[0])
        self.config["a"] = int(params[1])
        self.config["r_min"] = params[2]
        self.opt_counter += 1
        print "Iteration %s" % str(self.opt_counter)
        print "Model: K = %s,  a = %s  r_min = %.6E" % (str(self.config["k"]),
                                                        str(self.config["a"]),
                                                        self.config["r_min"])
        # Optimise bandwidths
        h_i, d_i = self.optimise_bandwidths()
        # Run smoothing
        self.run_smoothing(self.config["r_min"], h_i, d_i)
        print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
            np.sum(self.rates), np.min(self.rates), np.max(self.rates))
        # Analyse comparison with target catalogue
        probs = GridProbabilities(self.target_catalogue, self)
        # Count observed earthquakes
        probs.count_events()
        # Get Poisson Log-Likelihood
        poiss_llh = probs.poisson_loglikelihood()
        kagan_i0 = probs.get_i0()
        kagan_i1 = probs.get_i1()
        print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
            poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)
        return -poiss_llh#, kagan_i0, kagan_i1

class HelmstetterEtAl2007(HelmstetterNearestNeighbour):
    """
    Minimal adaptation of the Helmstetter methodology to consider only
    spatial smoothing
    """
    
    def build_distance_arrays(self):
        """
        
        """
        neq = self.learning_catalogue.get_number_events()
        distance_dset = self.fle.create_dataset("distance",
                                                (neq, neq),
                                                dtype="f")
        nearest_dset = self.fle.create_dataset("d_nearest",
                                               (neq, neq),
                                               dtype="i")
        counter1 = ProgressCounter(neq, 10.0, timer=True)
        for i in range(neq - 1):
            counter1.update(i)
            epi_distance = geodetic.geodetic_distance(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                self.learning_catalogue.data["longitude"],
                self.learning_catalogue.data["latitude"])
            distance_dset[i, :] = epi_distance
            nearest_dset[i, :] = np.argsort(epi_distance)
    
    def build_catalogue_2_grid_array(self):
        """
        In the hdf5 file store, for each earthquake, a Ngd by 4 array of
        dx1, dx2, dy1, dy2
        """
        neq = self.learning_catalogue.get_number_events()
        print "Building event to grid distance arrays"
        print str(datetime.now())
        gdist_dset = self.fle.create_dataset("grid_distance",
            (self.ngpts, 4, neq), chunks=(self.ngpts, 4, 1), dtype="f")
        print "Processing ..."
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        for i in range(neq):
            counter1.update(i)
            # Get distances
            dx1, dx2, dy1, dy2 = _get_distances(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                self.grid[:, 0], 
                self.grid[:, 1],
                self.grid_limits["xspc"],
                self.grid_limits["yspc"])
            gdist_dset[:, :, i] = np.column_stack([dx1, dx2, dy1, dy2])
        counter1.reset()
        print "done!"

    def optimise_bandwidths(self):
        """
        Nearest neighbour optimisation
        """
        # Build data from file
        neq = self.learning_catalogue.get_number_events()
        d_i = np.zeros(neq)
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        for i in range(neq):
            counter1.update(i)
            # Find Kth nearest neighbour
            locs = self.fle["d_nearest"][i, 1:(self.config["k"] + 1)]
            distances = self.fle["distance"][i, :]
            d_i[i] = np.max(distances[locs])
        d_i[d_i < 0.5] = 0.5
        return d_i
    
    def run_smoothing(self, r_min, d_i):
        """
        Runs the smoothing algorithm
        """
        self.rates = r_min * np.ones(self.ngpts)
        neq = self.learning_catalogue.get_number_events()
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        dtime = float(self.learning_catalogue.end_year -
                      self.learning_catalogue.start_year + 1)
        for i in range(self.learning_catalogue.get_number_events()):
            counter1.update(i)
            weight = self.learning_weights[i] / dtime
            dgeo = (self.fle["grid_distance"][:, :, i]).astype("float64")
            # Distance kernel
            # Get length scales based on bandwidth
            xls, yls = _get_length_scales(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                d_i[i])
            dx1 = dgeo[:, 0]
            dx2 = dgeo[:, 1]
            dy1 = dgeo[:, 2]
            dy2 = dgeo[:, 3]
            valid_long = np.logical_and(dx1 <= (5.92 * xls),
                                        dx2 >= (-5.92 * xls))
            valid_lat = np.logical_and(dy1 <= (5.92 * yls),
                                       dy2 >= (-5.92 * yls))
            idx = np.logical_and(valid_long, valid_lat)
            #print np.sum(valid_long), np.sum(valid_lat), np.sum(idx)
            distance_kernel = 0.25  * (
                erf(dx2[idx] / xls) - erf(dx1[idx] / xls)) *\
                (erf(dy2[idx] / yls) - erf(dy1[idx] / yls))
            self.rates[idx] += (weight * distance_kernel)
    
    def optimum_smoothing(self, params):
        """
        Runs the Helmstetter Adaptive Kernel Smoothing method optimized 
        """
        self.config["k"] = int(params[0])
        self.config["r_min"] = params[2]
        self.opt_counter += 1
        print "Iteration %s" % str(self.opt_counter)
        print "Model: K = %s,  r_min = %.6E" % (str(self.config["k"]),
                                                    self.config["r_min"])
        # Optimise bandwidths
        d_i = self.optimise_bandwidths()
        # Run smoothing
        self.run_smoothing(self.config["r_min"], d_i)
        print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
            np.sum(self.rates), np.min(self.rates), np.max(self.rates))
        # Analyse comparison with target catalogue
        probs = GridProbabilities(self.target_catalogue, self,
                                  self.target_weights)
        # Count observed earthquakes
        probs.count_events()
        # Get Poisson Log-Likelihood
        poiss_llh = probs.poisson_loglikelihood()
        uniform_llh = probs.uniform_spatial_poisson_loglikelihood()
        kagan_i0 = probs.get_i0()
        kagan_i1 = probs.get_i1()
        # Calculate probability gains last, as we mess with the catalogue
        if type(self.config["target_mmins"]) == list:
            prob_gain = probs.probability_gain()#self.config["target_mmins"],
                                               #self.config["bvalue"])
        else:
            self.config["target_mmins"] = None
            prob_gain = probs.probability_gain()#self.config["target_mmins"],
                                               #self.config["bvalue"])
#        print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f " \
#            "Uniform LLH = %.6f, Probability gain = %.6f"%(
#            poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1, uniform_llh, prob_gain)
        filename = 'error_diagram_b%.2f_%i_%i_%i_%i.png' %(
            self.config["bvalue"], self.learning_catalogue.start_year,
            self.learning_catalogue.end_year, self.target_catalogue.start_year,
            self.target_catalogue.end_year)
        probs.error_diagram(title='test',filename=filename)
        return -poiss_llh, kagan_i0, kagan_i1, uniform_llh, prob_gain

    def exhaustive_smoothing(self, krange, r_min_range):
        """

        """
        max_llh_params = [None, None]
        max_poiss_llh = -999999999999999.
        for kval in krange:
            for rval in r_min_range:
                self.config["k"] = kval
                self.config["r_min"] = rval
                print "Model: K = %s,  r_min = %.6E" % (str(self.config["k"]),
                                                        self.config["r_min"])
                # Optimise bandwidths
                d_i = self.optimise_bandwidths()
                # Run smoothing
                self.run_smoothing(self.config["r_min"], d_i)
                print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
                    np.sum(self.rates), np.min(self.rates), np.max(self.rates))
                # Analyse comparison with target catalogue
                probs = GridProbabilities(self.target_catalogue, self,
                                          self.target_weights)
                # Count observed earthquakes
                probs.count_events()
                # Get Poisson Log-Likelihood
                poiss_llh = probs.poisson_loglikelihood()
                kagan_i0 = probs.get_i0()
                kagan_i1 = probs.get_i1()
                print "Poisson LLH = %.6f" % poiss_llh #,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
#                    poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)
                if poiss_llh > max_poiss_llh:
                    max_llh_params = [kval, rval]
                    max_poiss_llh = poiss_llh
        return max_llh_params, max_poiss_llh

class FixedWidthSmoothing(HelmstetterEtAl2007):
    """
    Implements basic smoothing with a fixed with
    """
    def run_smoothing(self, r_min, d_i):
        """
        Runs the smoothing algorithm
        :param float d_i:
            Bandwidth km
        """
        self.rates = r_min * np.ones(self.ngpts)
        neq = self.learning_catalogue.get_number_events()
        counter1 = ProgressCounter(neq, 5.0, timer=True)
        dtime = float(self.learning_catalogue.end_year -
                      self.learning_catalogue.start_year + 1)
        for i in range(self.learning_catalogue.get_number_events()):
            counter1.update(i)
            weight = self.learning_weights[i] / dtime
            dgeo = (self.fle["grid_distance"][:, :, i]).astype("float64")
            # Distance kernel
            # Get length scales based on bandwidth
            xls, yls = _get_length_scales(
                self.learning_catalogue.data["longitude"][i],
                self.learning_catalogue.data["latitude"][i],
                d_i)
            dx1 = dgeo[:, 0]
            dx2 = dgeo[:, 1]
            dy1 = dgeo[:, 2]
            dy2 = dgeo[:, 3]
            valid_long = np.logical_and(dx1 <= (5.92 * xls),
                                        dx2 >= (-5.92 * xls))
            valid_lat = np.logical_and(dy1 <= (5.92 * yls),
                                       dy2 >= (-5.92 * yls))
            idx = np.logical_and(valid_long, valid_lat)
            #print np.sum(valid_long), np.sum(valid_lat), np.sum(idx)
            distance_kernel = 0.25  * (
                erf(dx2[idx] / xls) - erf(dx1[idx] / xls)) *\
                (erf(dy2[idx] / yls) - erf(dy1[idx] / yls))
            self.rates[idx] += (weight * distance_kernel)

    def optimum_smoothing(self, params):
        """
        Runs the Helmstetter Adaptive Kernel Smoothing method optimized 
        """
        self.config["bandwidth"] = int(params[0])
        self.config["r_min"] = params[2]
        self.opt_counter += 1
        print "Iteration %s" % str(self.opt_counter)
        print "Model: K = %s,  r_min = %.6E" % (str(self.config["bandwidth"]),
                                                    self.config["r_min"])
        # Run smoothing
        self.run_smoothing(self.config["r_min"], self.config["bandwidth"])
        print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
            np.sum(self.rates), np.min(self.rates), np.max(self.rates))
        # Analyse comparison with target catalogue
        probs = GridProbabilities(self.target_catalogue, self,
                                  self.target_weights)
        # Count observed earthquakes
        probs.count_events()
        # Get Poisson Log-Likelihood
        poiss_llh = probs.poisson_loglikelihood()
        uniform_llh = probs.uniform_spatial_poisson_loglikelihood()
        prob_gain = probs.probability_gain()
        kagan_i0 = probs.get_i0()
        kagan_i1 = probs.get_i1()
        print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
            poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)
        return -poiss_llh

    def exhaustive_smoothing(self, krange, r_min_range):
        """

        """
        for kval in krange:
            for rval in r_min_range:
                self.config["bandwidth"] = kval
                self.config["r_min"] = rval
                print "Model: Bandwitdh = %s km,  r_min = %.6E" % (
                    str(self.config["bandwidth"]), self.config["r_min"])
                # Run smoothing
                self.run_smoothing(self.config["r_min"],
                                   self.config["bandwidth"])
                print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
                    np.sum(self.rates), np.min(self.rates), np.max(self.rates))
                # Analyse comparison with target catalogue
                probs = GridProbabilities(self.target_catalogue, self,
                                          self.target_weights)
                # Count observed earthquakes
                probs.count_events()
                # Get Poisson Log-Likelihood
                poiss_llh = probs.poisson_loglikelihood()
                kagan_i0 = probs.get_i0()
                kagan_i1 = probs.get_i1()
                print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
                    poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)



class GridProbabilities(object):
    """
    Class to implement a suite of likelihood tests for comparing forecast
    and observed rates on a grid
    """
    def __init__(self, catalogue, grid, weights=None, mbins=None):
        """

        """
        self.catalogue = catalogue
        if catalogue.end_year:
            self.end_year = catalogue.end_year
        else:
            self.end_year = catalogue.update_end_year()

        if catalogue.start_year:
            self.start_year = catalogue.start_year
        else:
            self.start_year = catalogue.update_start_year()

        self.nyr = float(self.end_year - self.start_year + 1)

        self.grid = grid
        self.g_x, self.g_y = np.meshgrid(grid.lon_bins, grid.lat_bins)
        self.nlat, self.nlon = self.g_x.shape
        self.rates = np.reshape(np.copy(self.grid.rates),
                                [self.nlat - 1, self.nlon - 1])
        if mbins is not None:
            self.mbins = mbins
            self.n_m = len(mbins) - 1
        else:
            self.mbins = np.array([-np.inf, np.inf])
            self.n_m = 1
        self.obs_rates = np.zeros_like(self.rates)
        self.weights = weights
        self._get_area()
        
    def _subset_catalogue(self, mmin):
        """
        Given the configuration subset the catalogue
        """
       
        idx_m = np.greater_equal(
            self.catalogue.data["magnitude"],mmin)
        #print self.catalogue(
        self.catalogue.purge_catalogue(idx_m)
        self.weights = self.weights[idx_m] # fix here!!
#        self.catalogue = self.config["target_end"]
        
    def _get_area(self):
        """
        Gets the normalised cell area
        """
        self.area = np.zeros_like(self.rates)
        for i in range(self.nlat - 1):
            dlat = np.sin(np.radians(self.g_y[i + 1, 0])) -\
                np.sin(np.radians(self.g_y[i, 0]))
            self.area[i, :] = (np.pi / 180) * (6371.0 ** 2.) * dlat *\
                (self.g_x[0, 1:] - self.g_x[0, :-1])


    def count_events(self):
        """
        Counts the events in a given set of longitude, latitude and magnitude
        bins using historgramdd
        """
        xym = np.column_stack([self.catalogue.data["longitude"],
                               self.catalogue.data["latitude"],
                               self.catalogue.data["magnitude"]])
        if self.n_m == 1:
            self.obs_rates, _ = np.histogramdd(xym[:, :2],
                                               bins=[self.grid.lon_bins,
                                                     self.grid.lat_bins],
                                               weights=self.weights
                                               )
            self.obs_rates = self.obs_rates.T
        else:
            self.obs_rates, _ = np.histogramdd(
                xym, bins=[self.grid.lon_bins,
                           self.grid.lat_bins,
                           self.mbins],
                weights=self.weights)
        self.obs_rates /= float(self.catalogue.end_year -
                                self.catalogue.start_year + 1)

    def uniform_rate_density(self, nevents, nyrs):
        """
        Define an annual uniform rate density
        """
        annual_rate = float(nevents) / float(nyrs)
        # Get area for each cell
        return annual_rate * self.area / np.sum(self.area)

    def uniform_spatial_poisson_loglikelihood(self, mmin=0):
        """Returns loglikelood function for uniform 
        spatial distribution based on number of
        events in learning catalogue
        """
        
        idx_m = np.greater_equal(
            self.grid.learning_catalogue.data["magnitude"],mmin)
        self.grid.learning_catalogue.purge_catalogue(idx_m)
        nevents = self.grid.learning_catalogue.get_number_events()
        nyrs = float(self.grid.learning_catalogue.end_year - \
                         self.grid.learning_catalogue.start_year + 1)
        uniform_rates = self.uniform_rate_density(nevents,nyrs)
        # now calculate log likelihood of model
        #log_llh = np.sum(self.obs_rates*np.log(uniform_rates))
        log_llh = np.sum(np.log(((uniform_rates**self.obs_rates)* \
                                     np.exp(-uniform_rates))/factorial(self.obs_rates)))
        return log_llh

    def poisson_loglikelihood(self, rates=None):
        """
        Returns Poisson loglikelihood function
        """
        if rates is None:
            rates = self.rates
        idx = rates > 0.0
        return np.sum(np.log(
            (rates[idx] ** self.obs_rates[idx]) *\
                np.exp(-rates[idx]) / factorial(self.obs_rates[idx])))

    def probability_gain(self):
        """Probability gain as defined in equation 9 
        of Helmstetter et al. 2007.
        """
        self.avals = np.log10(self.rates) + self.grid.config['bvalue']*self.grid.config["mmin"]
        target_mmins = self.grid.config["target_mmins"]
        if target_mmins == None:
            poiss_llh = self.poisson_loglikelihood()
            uniform_llh = self.uniform_spatial_poisson_loglikelihood()
            print poiss_llh
            print uniform_llh
            prob_gain = np.exp((poiss_llh - uniform_llh)/ \
                                   self.catalogue.get_number_events())
        else:
            prob_gain = {}
            # sort to make purging the catalogue more straightforward
            target_mmins = sorted(target_mmins)
            for mmin in target_mmins:
                # Remove events < mmin from target and learning catalogue
                self._subset_catalogue(mmin)
                self.count_events() # Update target rates based on subsetted catalogue
                # Calculate gridded rates for M >= Mmin
                self.new_rates = np.power(10, 
                                          (self.avals - (self.grid.config["bvalue"]*mmin)))
                poiss_llh = self.poisson_loglikelihood(self.new_rates)
                uniform_llh = self.uniform_spatial_poisson_loglikelihood(mmin)
                print 'Mmin', mmin
                print 'poiss_llh', poiss_llh
                print 'uniform_llh', uniform_llh
                print 'number of events', self.catalogue.get_number_events()
                probability_gain = np.exp((poiss_llh - uniform_llh)/ \
                                              self.catalogue.get_number_events())
                print 'probability_gain', probability_gain
                prob_gain[str(mmin)]=probability_gain
        return prob_gain

    def get_i1(self):
        """
        Returns the Kagan (2009) I1 metric
        """

        # Total rate
        total_rate = float(np.sum(self.obs_rates))
        # Rate per km^2
        idx = self.obs_rates > 0.
        forecast_rate = self.rates[idx] / self.area[idx]
        rate_per_km2 = total_rate / np.sum(self.area)
        print 'total rate', total_rate
        print 'rate_per_km2', rate_per_km2
        print 'forecast rate', forecast_rate
        return (1. / total_rate) * np.sum(np.log2(forecast_rate /
                                                  rate_per_km2))

    def get_i0(self):
        """
        Returns the Kagan (2009) I0 metric
        """
        area_norm = self.area / np.sum(self.area)
        norm_rates = self.rates / np.sum(self.rates)
        idx = norm_rates > 0.0
        return np.sum(norm_rates[idx] * np.log2(norm_rates[idx] /
                                                area_norm[idx]))

    def error_diagram(self, title=None, filename=None, filetype="png", dpi=300):
        """
        Creates the Molchan area diagram
        """
        # Normalize forecast and observed rates
        forecast = (self.rates / np.sum(self.rates)).flatten()
        observed = (self.obs_rates / np.sum(self.obs_rates)).flatten()
        area_norm = (self.area / np.sum(self.area)).flatten()
        #area_norm = (area / np.sum(area)).flatten()
        num_cell = len(forecast)
        
        # Poisson Model
        poisson_idx = np.argsort(area_norm)[::-1]
        poisson = 1.0 - np.cumsum(area_norm[poisson_idx])
        area_poisson = np.cumsum(area_norm[poisson_idx])
        
        # Forecase rate
        forecast_idx = np.argsort(forecast)[::-1]
        cum_forecast = 1.0 - np.cumsum(forecast[forecast_idx])
        area_forecast = np.cumsum(area_norm[forecast_idx])
        
        # Observed Rate
        observed_idx = np.argsort(observed)[::-1]
        cum_observed = 1.0 - np.cumsum(observed[observed_idx])
        area_observed = np.cumsum(area_norm[observed_idx])
        # Plot
        plt.figure(figsize=(8,8))
        plt.semilogx(area_poisson, poisson, 'k-', lw=2, label="Poisson")
        plt.semilogx(area_forecast, cum_forecast, 'r-', lw=2, label="Model")
        plt.semilogx(area_observed, cum_observed, 'b-', lw=2, label="Observed")
        plt.ylim(0., 1.)
        plt.yticks(np.arange(0., 1.1, 0.1))
        plt.xlabel(r"Fraction of Alarm Area, $\tau$", fontsize=16)
        plt.ylabel(r"Fraction of Failures to Predict, $\nu$", fontsize=16)
        plt.legend(loc=3, fontsize=16, framealpha=1.)
        if title:
            plt.title(title, fontsize=18)
        plt.grid(True)
        if filename:
            plt.savefig(filename, format=filetype, dpi=dpi)


def optimum_smoothing(params, model):
    """
    Runs the Helmstetter Adaptive Kernel Smoothing method optimized 
    """
    model.config["k"] = params[0]
    model.config["a"] = params[1]
    model.config["r_min"] = params[2]
    print "Model: K = %s,  a = %s  r_min = %.6E" % (str(model.config["k"]),
                                                    str(model.config["a"]),
                                                    model.config["r_min"])
    # Optimise bandwidths
    h_i, d_i = model.optimise_bandwidths()
    # Run smoothing
    model.run_smoothing(model.config["r_min"], h_i, d_i)
    print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (np.sum(model.rates),
                                                           np.min(model.rates),
                                                           np.max(model.rates))
    # Analyse comparison with target catalogue
    probs = GridProbabilities(model.target_catalogue, model)
    # Count observed earthquakes
    probs.count_events()
    # Get Poisson Log-Likelihood
    poiss_llh = probs.poisson_loglikelihood()
    kagan_i0 = probs.get_i0()
    kagan_i1 = probs.get_i1()
    print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
        poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)
    return -poiss_llh

def run_helmstetter_full(catalogue, bbox, config, completeness,
        output_file_stem, tmp_file="tmp1.hdf5"):
    """
    Executes the full function
    """
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    model1 = HelmstetterNearestNeighbour(bbox, config, catalogue,
                                         storage_file=tmp_file)
    # Set up catalogues
    model1._get_catalogue_completeness_weights(completeness)
    print model1.learning_catalogue.get_number_events()
    # Build catalogue arrayree
    model1.build_time_distance_arrays()
    # Build grid array
    model1.build_catalogue_2_grid_array(float(config["target_start"]),
                                        float(config["ndays"]))

    # Solver with OPENOPT
    initial_params = [5, 300, 1.0E-4]
    solver = GLP(model1.optimum_smoothing, initial_params, maxiter=1E3)
    solver.lb = [1, 0, 0.0]
    solver.ub = [100, 1000, 1.0]
    solver.name = "glp_1"
    solver.discrtol = 1.1E-5
    results = solver.minimize("de", plot=False)
    # Now run optimisation using L-BFGS-B
    # Bounds
    #bnds = [(0.0, None), (0.0, None), (0.0, None)]

    #results = minimize(optimum_smoothing,
    #                   initial_params,
    #                   args=[model1],
    #                   method="Nelder-Mead")
    #if not results.success:
    #    print results.message
    #    raise ValueError("Algorithm Failed :'(!")

    config["k"], config["a"], config["r_min"] = results.xf
    print "Best Results"
    print "K = %s  a = %s  rmin = %.8e" % (str(config["k"]),
                                           str(config["a"]),
                                           config["r_min"])
    # Run algorithm with optimum parameters
    h_i, d_i = model1.optimise_bandwidths()
    model1.run_smoothing(config["r_min"], h_i, d_i)
    # Export grid and rates
    outputs = np.column_stack([model1.grid, model1.rates])
    np.savetxt(output_file_stem + ".csv", outputs, fmt="%.3f,%.3f,%.8e")
    # Plot and Save Molhan Diagram
    probs = GridProbabilities(model1.target_catalogue, model1)
    probs.count_events()
    error_diag_name = output_file_stem + "Error_Diagram.png"
    probs.error_diagram(filename=error_diag_name,
                        filetype="png", dpi=300)


def run_helmstetter_spatial(catalogue, bbox, config, completeness,
        output_file_stem, comp_file, tmp_file="tmp1.hdf5"):
    """
    Executes the full function
    """
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    model1 = HelmstetterEtAl2007(bbox, config, catalogue,
                                         storage_file=tmp_file)
    # Set up catalogues
    #model1._get_catalogue_completeness_weights(completeness)
    model1._get_spatio_temporal_completeness_weights(comp_file)
    print model1.learning_catalogue.get_number_events()
    # Build catalogue arrayree
    model1.build_distance_arrays()
    # Build grid array
    model1.build_catalogue_2_grid_array()
    # Optimise bandwitdhs
    d_i = model1.optimise_bandwidths()
    model1.run_smoothing(config["r_min"], d_i)
    print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
        np.sum(model1.rates), np.min(model1.rates), np.max(model1.rates))

    # Export grid and rates
    outputs = np.column_stack([model1.grid, model1.rates])
    np.savetxt(output_file_stem + ".csv", outputs, fmt="%.3f,%.3f,%.8e")
    # Plot and Save Molchan Diagram
    probs = GridProbabilities(model1.target_catalogue, model1)
    probs.count_events()
    poiss_llh = probs.poisson_loglikelihood()
    kagan_i0 = probs.get_i0()
    kagan_i1 = probs.get_i1()
    print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
        poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)

    #error_diag_name = output_file_stem + "Error_Diagram.png"
    #probs.error_diagram(filename=error_diag_name,
    #                    filetype="png", dpi=300)

def run_fixed_width_spatial(catalogue, bbox, config, completeness,
        output_file_stem, comp_file, tmp_file="tmp1.hdf5"):
    """
    Executes the full function
    """
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    model1 = FixedWidthSmoothing(bbox, config, catalogue,
                                 storage_file=tmp_file)
    # Set up catalogues
    #model1._get_catalogue_completeness_weights(completeness)
    model1._get_spatio_temporal_completeness_weights(comp_file)
    print model1.learning_catalogue.get_number_events()
    # Build catalogue arrayree
    model1.build_distance_arrays()
    # Build grid array
    model1.build_catalogue_2_grid_array()
    # Optimise bandwitdhs
    d_i = model1.optimise_bandwidths()
    model1.run_smoothing(config["r_min"], d_i)
    print "Total Rates: %.6f  (min = %.8f, max = %.8f)" % (
        np.sum(model1.rates), np.min(model1.rates), np.max(model1.rates))

    # Export grid and rates
    outputs = np.column_stack([model1.grid, model1.rates])
    np.savetxt(output_file_stem + ".csv", outputs, fmt="%.3f,%.3f,%.8e")
    # Plot and Save Molchan Diagram
    probs = GridProbabilities(model1.target_catalogue, model1)
    probs.count_events()
    poiss_llh = probs.poisson_loglikelihood()
    kagan_i0 = probs.get_i0()
    kagan_i1 = probs.get_i1()
    print "Poisson LLH = %.6f,  I0 = %.6f,   I1 = %.6f,   I' = %.6f" %(
        poiss_llh, kagan_i0, kagan_i1, kagan_i0 - kagan_i1)






SARA_COMP_TABLE = np.array([[1992., 4.5],
                            [1974., 5.],
                            [1964., 5.5],
                            [1954., 5.75],
                            [1949., 6.],
                            [1949., 6.5],
                            [1930., 7.0]])
SARA_CAT_FILE = "catalogue/sara_all_v07_harm_per123_dist_crustal_clean.csv"
SARA_DECLUST_CAT = "catalogue/sara_cat_shallow_declust.csv"
COMP_FILE = "sam_completeness_zones.hdf5"


if __name__ == "__main__":
    # Load in catalogue
    parser = CsvCatalogueParser(SARA_DECLUST_CAT)
    cat1 = parser.read_file()
    idx = cat1.data["magnitude"] >= 3.0
    cat1.purge_catalogue(idx)
    bbox = [-90.5, -30.0, 0.1, -60.5, 15.5, 0.1, 0., 100.0, 100.0]
    config = {"k": 5,
              "r_min": 0.0005 / 25.0, 
              "bvalue": 1.0, "mmin": 3.0,
              "learning_start": 1930, "learning_end": 2003,
              "target_start": 2004, "target_end": 2013}
    # Run
    run_helmstetter_spatial(cat1, bbox, config, SARA_COMP_TABLE, 
                            "helmstetter_SA_trial_v3",
                            COMP_FILE,
                            tmp_file="hsa_tmp.hdf5")


                
