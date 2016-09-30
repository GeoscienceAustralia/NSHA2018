"""Script to convert between fault slip-rates and Gutenberg-Richter a and b values
Jonathan Griffin, AIFDR
October 2012
Modified Jonathan Griffin, June 2016
"""

import math
import  numpy as np

def slip2GR(slip_rate, A, b, M_max, M_min=0, mu=3e11):
    """Calculate Gutenberg-Richter a values
    from fault geometry and slip-rate based on Youngs
    and Coppersmith 1985
    :param slip_rate:
        Float: Fault slip-rate in mm per year
    :param A:
        Float: Fault area in km^2
    :param b:
        Float: Gutenberg-Richter b value
    :param M_max:
        Maximum magnitude for Gutenberg-Richter 
        distribution
    :param M_min:
        Float: Minimum magnitude for which
        'a' is calculated (normally 0)
    :param mu:
        Float: Crustal shear modulus in dyn/cm^2
    :returns a:
        Float: Gutenberg-Richter a value"""

    c = 1.5 # Hanks and Kanamori coefficient
    beta = np.log(10)*b
    slip_fault = slip_rate/10  # convert to cm/a
    A = A*np.power(10,10) # convert to cm^2
    # Convert moment magnitude to moment for M_max
    moment_max = np.power(10,(c*M_max + 16.05)) # dyn/cm = N/km
    
    # Calculate number of eqs > M_min
    N = (mu*A*(slip_fault)*(c-b)*(1-np.exp(-1*beta*(M_max-M_min))))/(b*moment_max*np.exp(-1*beta*(M_max-M_min)))
    # Calculate a
    a = np.log10(N) + b*M_min
    # Momement rate from slip rate
    moment_rate = mu*A*slip_fault #dyncm
    moment_rate = moment_rate/1e7 #Nm

    return a, moment_rate

def GR2sliprate(a, b, A, M_max, M_min=0, mu=3e11):
    """Calculate slip-rates for a fault from 
    Gutenberg-Richter a values based on Youngs
    and Coppersmith 1985
    :param a:
        Float: Gutenberg-Richter a value
    :param A:
        Float: Fault area in km^2
    :param b:
        Float: Gutenberg-Richter b value
    :param M_max:
        Maximum magnitude for Gutenberg-Richter 
        distribution
    :param M_min:
        Float: Minimum magnitude for which
        'a' is defined (normally 0)
    :param mu:
        Float: Crustal shear modulus in dyn/cm^2
    :returns slip_rate:
        Float: Fault slip-rate in mm per year
    :returns moment_rate
        Float: Seismic moment rate across fault"""

    c = 1.5 # Hanks and Kanamori coefficient
    beta = np.log(10)*b
    A = A*np.power(10,10) # convert to cm^2
    # Convert moment magnitude to moment for M_max
    moment_max = np.power(10,(c*M_max + 16.05)) # dyn/cm = N/km
    
    # Convert a to N, number of eqs > M_min
    N = np.power(10,a - b*M_min)
    # Calculate slip_rate (cm/a)
    slip_rate = (b*N*moment_max*np.exp(-1*beta*(M_max-M_min)))/(mu*A*(c-b)*(1-np.exp(-1*beta*(M_max-M_min))))
    slip_rate = slip_rate*10 # Convert to mm/a
    moment_rate = (b*N*moment_max*np.exp(-1*beta*(M_max-M_min)))/((c-b)*(1-np.exp(-1*beta*(M_max-M_min))))
    moment_rate = moment_rate/1e7 #Nm
    return slip_rate, moment_rate

if __name__ == "__main__":

    # Recurrence parameters calculated from fitting seimicity data to a G-R curve
    M_min = 4.9 # Minimum magnitude (for catalogue copmmleteness and G-R generation
    M_max = 7.0 # Maximum magnitude in source zone
    lambda_min =  1.34 #0.175 ##2.925 Number of earthquakes greater than M_min per year (= A_min)
    b = 0.94 # Gutenberg-Richter b-value
    beta = 2.303*b

    # Fault Parameter
    slip_fault = 30.0 # Reported slip rate for the fault (mm/a)
    L = 370 #56. #146 #93 #299    # fault length (km)
    W = 45. # 18. # 14.8 #18. .#18 # fault width
    A = L*W # Fault area




