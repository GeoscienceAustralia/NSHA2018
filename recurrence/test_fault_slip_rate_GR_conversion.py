"""Unit test for conversion of fault slip rate to 
Gutenberg-Richter a values. Based on formulation of
Youngs and Coppersmith 1985
"""

import unittest
from numpy.testing import assert_allclose
from fault_slip_rate_GR_conversion import slip2GR, GR2sliprate

class slip2GRTestCase(unittest.TestCase):
    
    def setUp(self):
        self.test_sliprates = [0.1, 0.1, 3., 3., 3., 3.]
        self.test_areas = [200., 2000., 2000., 2000., 2000., 2000.]
        self.test_bvalues = [1.0, 1.0, 1.0, 1.0, 0.9, 1.1]
        self.test_mmaxs = [6.5, 6.5, 6.5, 6.5, 6.5, 7.2]
        self.test_mmins = [0, 0, 0, 4., 0., 0.]
        self.test_mus = [3e11, 3e11, 3e11, 3e11, 3e11, 3e11]

        self.expected_avalue = [2.17712, 3.17712, 4.65424, 4.652867, 4.12918, 4.885939]
        self.expected_moment_rate = [6e14, 6e15, 1.8e17, 1.8e17, 1.8e17, 1.8e17]

    def test_avalues(self):
        for i in range(len(self.test_sliprates)):
            avalue, moment_rate = slip2GR(self.test_sliprates[i],
                                          self.test_areas[i],
                                          self.test_bvalues[i],
                                          self.test_mmaxs[i],
                                          self.test_mmins[i],
                                          self.test_mus[i])
            self.assertAlmostEqual(self.expected_avalue[i],
                                   avalue, places = 5)

    def test_momentrates(self):
        for i in range(len(self.test_sliprates)):
            avalue, moment_rate = slip2GR(self.test_sliprates[i],
                                          self.test_areas[i],
                                          self.test_bvalues[i],
                                          self.test_mmaxs[i],
                                          self.test_mmins[i],
                                          self.test_mus[i])
            self.assertAlmostEqual(self.expected_moment_rate[i],
                                   moment_rate, places = 5)

class GR2sliprateTestCase(unittest.TestCase):
    def setUp(self):
        self.test_areas = [200., 2000., 2000., 2000., 2000., 2000.]
        self.test_bvalues = [1.0, 1.0, 1.0, 1.0, 0.9, 1.1]
        self.test_mmaxs = [6.5, 6.5, 6.5, 6.5, 6.5, 7.2]
        self.test_mmins = [0, 0, 0, 4., 0., 0.]
        self.test_mus = [3e11, 3e11, 3e11, 3e11, 3e11, 3e11]
        self.test_avalues = [2.1771211174, 3.1771211174, 4.6542423721, 4.6528669736,
                             4.1291806326, 4.8859398061]

        self.expected_sliprates = [0.1, 0.1, 3., 3., 3., 3.]
        self.expected_moment_rate = [6e14, 6e15, 1.8e17, 1.8e17, 1.8e17, 1.8e17]

    def test_sliprates(self):
        for i in range(len(self.test_avalues)):
            sliprate, moment_rate = GR2sliprate(self.test_avalues[i],
                                                self.test_bvalues[i],
                                                self.test_areas[i],
                                                self.test_mmaxs[i],
                                                self.test_mmins[i],
                                                self.test_mus[i])
            self.assertAlmostEqual(self.expected_sliprates[i],
                                   sliprate, places = 5)

    def test_momentrates(self):
        for i in range(len(self.test_avalues)):
            avalue, moment_rate = GR2sliprate(self.test_avalues[i],
                                          self.test_bvalues[i],
                                          self.test_areas[i],
                                          self.test_mmaxs[i],
                                          self.test_mmins[i],
                                          self.test_mus[i])
            assert_allclose(self.expected_moment_rate[i],
                                   moment_rate, rtol = 1e-5)
