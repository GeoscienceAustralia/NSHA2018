# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2013-2017 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Module exports :class:`GaullEtAL1990PGA`,
"""
from __future__ import division

import numpy as np
from scipy.constants import g

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA

"""
Do Ghasemi 2017 fixed-hinge bilinear GOR ML2MW conversion (in reverse) as 
used for the Australian 2018 National Seismic Hazard Assessment as referred to in:

Allen, T., J. Griffin, D. Clark, H. Ghasemi, L. Leonard, and T. Volti (2017). 
Towards the 2018 National Seismic Hazard Assessment: Draft design values as 
proposed for the Standards Australia AS1170.4–2018 Earthquake Design Actions, 
Geoscience Australia Professional Opinion 2017/02, Canberra, 44 pp.
"""
def ghasemi_bl_mw2ml(mw):
    a1 = 0.661
    a2 = 1.209
    a3 = 0.987
    hx = 4.25
    hy = a1 * hx + a2
    
    ml = (mw - a2) / a1
    
    if ml > hx:
        ml = hx + (mw - hy) / a3
        
    return ml

class GaullEtAL1990PGA(GMPE):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    """
    #: Supported tectonic region type is stable continental crust
    #: that the equations have been derived for Australia
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
    ])

    #: Supported intensity measure component is vertical
    #: :attr:`~openquake.hazardlib.const.IMC.VERTICAL`,
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.VERTICAL

    #: Supported standard deviation type is total
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: site params are not required
    REQUIRES_SITES_PARAMETERS = set()

    #: Required rupture parameter is magnitude
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', ))

    #: Required distance measure is Rjb distance
    #: see paragraph 'Predictor Variables', page 6.
    REQUIRES_DISTANCES = set(('rhypo', ))
    
    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        C = self.COEFFS[imt]
        
        magML = ghasemi_bl_mw2ml(rup.mag)
        c0 = 0.0
        R = np.sqrt(dists.rhypo ** 2 + c0 ** 2)
                
        mean = C['a'] * np.exp(C['b'] * magML) * R**(-1 * C['c'])
        
        if isinstance(imt, (PGA)):
            # For PGA convert from m/s**2 to ln g
            mean = np.log(mean / g)
        else:
            # For PGV convert from mm/s to cm/s/s
            mean = np.log(mean / 10.)

        stddevs = self._get_stddevs(C, stddev_types,  dists.rhypo.shape[0])
        
        return mean, stddevs

    def _get_stddevs(self, C, stddev_types, num_sites):
        """
        Return total standard deviation.
        """
        assert all(stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                   for stddev_type in stddev_types)
        stddevs = [np.zeros(num_sites) + C['sigma'] for _ in stddev_types]
        return stddevs

    #: dummy coefficient table 
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a     b      c      sigma
    pga  1     1      1      1
   
    """)


class GaullEtAL1990WesternAustralia(GaullEtAL1990PGA):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Coefficients for Western Australia
    """
    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  0.025  1.10   1.03   0.28
    pgv  3.30   1.04   0.96   0.28
    """)

class GaullEtAL1990SoutheasternAustralia(GaullEtAL1990PGA):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Coefficients for Southeastern Australia
    """

    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  0.088  1.10   1.20   0.28
    pgv  12.20  1.04   1.18   0.28
    """)
    
class GaullEtAL1990NortheasternAustralia(GaullEtAL1990PGA):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Coefficients for Northeastern Australia
    """
    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  0.06  1.04   1.08   0.32
    pgv  7.28  0.97   1.01   0.32
    """)
    
class GaullEtAL1990Indonesia(GaullEtAL1990PGA):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Coefficients for Indonesia
    """
    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  37.6   1.11   2.075   0.16
    pgv  3019.  1.04   1.940   0.16
    """)


class GaullEtAL1990PGAfromPGV(GMPE):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Uses PGV coefficients and divides velocity (in mm/s) by 750 to obtain 
    PGA in g as per Applied Technology Council (1984). Tentative provisions for 
    the development of seismic regulations for buildings, Applied Technology 
    Council ATC 3-06 Amended, Redwood City, CA, 505 pp.
    """
    #: Supported tectonic region type is stable continental crust
    #: that the equations have been derived for Western North America
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL

    #: Supported intensity measure types are spectral acceleration,
    #: and peak ground acceleration
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
    ])

    #: Supported intensity measure component is vertical
    #: :attr:`~openquake.hazardl'+str('%0.2f' % mw)ib.const.IMC.VERTICAL`,
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.VERTICAL

    #: Supported standard deviation type is total
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL
    ])

    #: site params are not required
    REQUIRES_SITES_PARAMETERS = set()

    #: Required rupture parameter is magnitude
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', ))

    #: Required distance measure is Rjb distance
    #: see paragraph 'Predictor Variables', page 6.
    REQUIRES_DISTANCES = set(('rhypo', ))

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See :meth:`superclass method
        <.base.GroundShakingIntensityModel.get_mean_and_stddevs>`
        for spec of input and result values.
        """
        C = self.COEFFS[imt]

        magML = ghasemi_bl_mw2ml(rup.mag)
        c0 = 0.0
        R = np.sqrt(dists.rhypo ** 2 + c0 ** 2)
        
        mean = C['a'] * np.exp(C['b'] * magML) * R**(-1 * C['c'])
        
        """
        # convert from mm/s to ln g using approximate ATC (1984) conversions
        # as assumed in: 
        
        McCue, K. (1993). The revised Australian seismic hazard map, 1991, 
        Proceedings of the 1993 Australian Earthquake Engineering Society 
        Conference, Melbourne, VIC
        """
        mean = np.log(mean / 750.)

        stddevs = self._get_stddevs(C, stddev_types,  dists.rhypo.shape[0])

        return mean, stddevs

    def _get_stddevs(self, C, stddev_types, num_sites):
        """
        Return total standard deviation.
        """
        assert all(stddev_type in self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
                   for stddev_type in stddev_types)
        stddevs = [np.zeros(num_sites) + C['sigma'] for _ in stddev_types]
        return stddevs

    #: dummy coefficient table 
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a     b      c      sigma
    pga  1     1      1      1
   
    """)


class GaullEtAL1990PGAfromPGVWesternAustralia(GaullEtAL1990PGAfromPGV):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Uses PGV coefficients and divides velocity (in mm/s) by 750 to obtain 
    PGA in g as per Applied Technology Council (1984). Tentative provisions for 
    the development of seismic regulations for buildings, Applied Technology 
    Council ATC 3-06 Amended, Redwood City, CA, 505 pp.
    
    Coefficients for Western Australia
    """

    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  3.30   1.04   0.96   0.28
    """)


class GaullEtAL1990PGAfromPGVSoutheasternAustralia(GaullEtAL1990PGAfromPGV):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Uses PGV coefficients and divides velocity (in mm/s) by 750 to obtain 
    PGA in g as per Applied Technology Council (1984). Tentative provisions for 
    the development of seismic regulations for buildings, Applied Technology 
    Council ATC 3-06 Amended, Redwood City, CA, 505 pp.
    
    Coefficients for Southeastern Australia
    """

    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  12.20  1.04   1.18   0.28
    """)
    
class GaullEtAL1990PGAfromPGVNortheasternAustralia(GaullEtAL1990PGAfromPGV):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Uses PGV coefficients and divides velocity (in mm/s) by 750 to obtain 
    PGA in g as per Applied Technology Council (1984). Tentative provisions for 
    the development of seismic regulations for buildings, Applied Technology 
    Council ATC 3-06 Amended, Redwood City, CA, 505 pp.
    
    Coefficients for Northeastern Australia
    """

    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  7.28  0.97   1.01   0.32
    """)
    
class GaullEtAL1990PGAfromPGVIndonesia(GaullEtAL1990PGAfromPGV):
    """
    Implement equation used to underpin the Structural design actions, part 4: 
    Earthquake actions in Australia, Standards Australia AS 1170.4–2007.
    Equation coefficients as per Gaull, B. A., M. O. Michael-Leiba, and 
    J. M. W. Rynn (1990). Probabilistic earthquake risk maps of Australia, 
    Aust. J. Earth. Sci. 37, 169-187, doi: 10.1080/08120099008727918.
    
    Uses PGV coefficients and divides velocity (in mm/s) by 750 to obtain 
    PGA in g as per Applied Technology Council (1984). Tentative provisions for 
    the development of seismic regulations for buildings, Applied Technology 
    Council ATC 3-06 Amended, Redwood City, CA, 505 pp.
    
    Coefficients for Indonesia
    """

    #: coefficient table as per Table 4 in Gaull et. al. (1990)
    COEFFS = CoeffsTable(sa_damping=5, table="""\
    IMT  a      b      c      sigma
    pga  3019.  1.04   1.940   0.16
    """)

