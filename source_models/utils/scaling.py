"""Implement other iterations of Leonard 2014 scaling

Reference: Leonard, M., 2014. Self-consistent earthquake fault-scaling relations:
Update and extension to stable continental strike-slip faults.
Bulletin of the Seismological Society of America, 104(6), pp 2953-2965.

Jonathan Griffin
Geoscience Australia May 2018
"""

from numpy import power, log10
from openquake.hazardlib.scalerel import Leonard2014_SCR

class Leonard2014_SCR_extension(Leonard2014_SCR):
    """
    Leonard, M., 2014. Self-consistent earthquake fault-scaling relations:
    Update and extension to stable continental strike-slip faults.
    Bulletin of the Seismological Society of America, 104(6), pp 2953-2965.

    Extends OpenQuake implementation to include
    length and width based scaling relationships (not just area and moment).
    """

    def get_width_from_length(self, length, rake):
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike-slip
            if length <=1600:
                return power(10,(1.0*log10(length)))
            elif length > 1600 and length <= 70000:
                return power(10,(1.068 + 0.667*log10(length)))
            else:
                return power(10,4.298)
        else:
            # Dip-slip (thrust or normal)
            if length <= 2500:
                return power(10, (1.0 + 1.0*log10(length)))
            else:
                return power(10, (1.130 + 0.667*log10(length)))

    def get_width_from_area(self, area, rake):
        if (-45 <= rake <= 45) or (rake >= 135) or (rake <= -135):
            # strike-slip
            if area <=2.56e6:
                return power(10, (0.5*log10(area)))
            elif area > 2.56e6 and area < 1400e6:
                return power(10,(0.641 + 0.4*log10(area)))
            else:
                return power(10,4.298)
        else:
            # Dip-slip (thrust or normal)
            if area <= 6.2e6:
                return power(10, (0.5*log10(area)))
            else:
                return power(10, (0.678 + 0.4*log10(area)))
