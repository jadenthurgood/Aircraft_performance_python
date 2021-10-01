import numpy as np
import os

def interp(h):
    """This function takes a height in ft and returns the geopotential height in ft,
    the temperature (R), Pressure (lbf/ft^2), Density (slugs/ft^3) and Speed of Sound (ft/s)

    Args:
        h (float): Altitude (MSL)

    Returns:
        tuple: (geopotential height, temperature, pressure, density, speed of sound)
    """
    
    #import the English Standard Atmosphere table of data from the CSV file
    #data comes from Phillips "Mechanics of Flight" Appendix A
    data = np.genfromtxt("./Standard_Atmosphere/Standard_Atmosphere_English_table.csv",dtype='float',delimiter=',',skip_header=1)
    
    geo = np.interp(h,data[:,0],data[:,1])
    temp = np.interp(h,data[:,0],data[:,2])
    pres = np.interp(h,data[:,0],data[:,3])
    den = np.interp(h,data[:,0],data[:,4])
    mach = np.interp(h,data[:,0],data[:,5])
    return (geo,temp,pres,den,mach)