"""Functions for calculating local magnitudes dependent on the 
location of the earthquake

Jonathan Griffin
Geoscience Australia, April 2016"""

from numpy import log10
from shapely.geometry import Polygon, Point

def ml_swa(amp, distance):
    """Calculate local magnitude for South Western Australia
    :param amp: 
        Float: Wood-Anderson amplitude in mm
    :param distance:
        Float: Distance from earthquake hypocentrer in km
        to station
    :returns: Local magnitude"""
    
    ml = log10(amp) + 1.137*log10(distance) + 0.000657*distance + 0.66
    return ml
    
def ml_sa(amp, distance):
    """Calculate local magnitude for South Australia
    :param amp: 
        Float: Wood-Anderson amplitude in mm
    :param distance:
        Float: Distance from earthquake hypocentrer in km
        to station
    :returns: Local magnitude"""
    
    ml = log10(amp) + 1.1*log10(distance) + 0.0013*distance + 0.7
    return ml    

def ml_sea(amp, distance):
    """Calculate local magnitude for South Eastern Australia
    :param amp: 
        Float: Wood-Anderson amplitude in mm
    :param distance:
        Float: Distance from earthquake hypocentrer in km
        to station
    :returns: Local magnitude"""
    
    ml = log10(amp) + 1.34*log10(distance/100.) + 0.00055*(distance-100) + 3.13
    return ml

def ml_richter(amp, distance):
    """calculate original Richter magnitude for earthquakes
    outside of Australia
    param amp: 
        Float: Wood-Anderson amplitude in mm
    :param distance:
        Float: Distance from earthquake hypocentrer in km
        to station
    :returns: Local magnitude"""
    pass

def calculate_local_magnitude(wa_amp, hypocentre, distance):
    """Calculate local magnitude dependent on location
    
    :param wa_amp:
        Float. Wood-Anderson zero to peak amplitude (mm)
    :param hypocentre:
        list: [longitude, latitude, depth(km)] of hypocentre
    :param distance:
        float: distance from hypocentre to station
    :returns: Local magnitude (Ml)"""
    
    # Define regions for differen ML calculations
    # Regions are sourced from http://rhe-eqm-swdev.dev.lan:8000/trac/wiki/project_Ml
    sw_australia_poly = Polygon(((110,-40), (110,-11), (140,-11), (140,-40)))
    south_australia_poly = Polygon(((130,-40), (130, -25), (140, -25), (140, -40)))
    se_australia_poly = Polygon(((140,-45), (140,-10), (155,-10), (155,-45)))
    # Heirarchicly test if point is in polygon and calculate relevant magnitude
    # Following ATWS calculations, defaults to standard Richter coefficients if 
    # not within Australian region
    eq_point = Point((hypocentre[0], hypocentre[1]))
    if eq_point.within(south_australia_poly):
        local_magnitude = ml_sa(wa_amp, distance)
    elif eq_point.within(sw_australia_poly):
        local_magnitude = ml_swa(wa_amp, distance)
    elif eq_point.within(se_australia_poly):
        local_magnitude = ml_sea(wa_amp, distance)
    else:
        pass
        #local_magnitude = ml_richter(wa_amp, distance)
        
    try:
        return local_magnitude
    except NameError:
        msg = 'Earthquake occurs outside of Australian region. ' \
            'No local magnitude scale defined'
        raise Exception, msg
        
    
