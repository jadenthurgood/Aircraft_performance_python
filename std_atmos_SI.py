import numpy as np

def interp(h):
    """ This function takes a height in meters and returns the geopotential height in meters,
    the temperature (K), Pressure (N/m^2), Density (kg/m^3) and Speed of Sound (m/s) """
   
    #import the SI Standard Atmosphere table of data from the CSV file
    data = np.genfromtxt("Standard_Atmosphere_SI_table.csv",dtype='float',delimiter=',',skip_header=1)
    
    geo = np.interp(h,data[:,0],data[:,1])
    temp = np.interp(h,data[:,0],data[:,2])
    pres = np.interp(h,data[:,0],data[:,3])
    den = np.interp(h,data[:,0],data[:,4])
    mach = np.interp(h,data[:,0],data[:,5])
    return (geo,temp,pres,den,mach)