import numpy as np
import std_atmos_Eng
import std_atmos_SI
import warnings

class Airplane():
    def __init__(self):
        self.wingspan = 10  #ft
        self.c_bar = 2 #ft
        self.T0 = 5000 #lbf
        self.AR = self.wingspan/self.c_bar
        self.e = 0.8 #Oswald efficiency
        self.Sw = self.wingspan*self.c_bar #wing area
        self.units = 'English'

    def get_density(self,alt,**kwargs):
        """ This function accepts the altitude as the argument and returns the density in slugs/ft^3 for "English" and N/m^3 for "SI"
        The user may specify the units in the input file for the class or using the keyword "units" and specifying "English" or "SI".
        If the user specifies the units using the keyword in the function call, this will overrite the units set in the input file. """

        if kwargs: #check if the user input kwargs
            if 'units' in kwargs: #check if the user used the kwarg 'units'
                if kwargs['units'] == 'SI': #check if the user specified the units at the time of function call. This will overrite the Input file
                    rho = std_atmos_SI.interp(alt)[3] #return the density from the interpolation
                elif kwargs['units'] == 'English':
                    rho = std_atmos_Eng.interp(alt)[3]
                else: #stop the code if the user has manually input an incorrect selection for units. 
                    raise ValueError("The only keyword argument accepted values for 'units' are 'SI' and 'English'.")
            else: #if the user did not use 'units' as a kwarg then raise value error
                raise ValueError("The only keyword argument used in get_density is 'units' which can be assigned the value 'SI' or 'English'.")
        
        else:  #check which units you are using in the input file
            if self.units == 'English':
                rho = std_atmos_Eng.interp(alt)[3]  
            else: #Use SI
                if self.units != 'SI': #Should be SI if not English. If not, then print a warning for the user letting them know defaults are SI, but still return density assuming SI
                    warnings.warn("Units were not specified in the input file as 'English' or 'SI' and no inputs were given in the function call. Defaulting to SI")
                rho = std_atmos_SI.interp(alt)[3]
        
        return rho

    def get_AR(self):
        return self.wingspan/self.c_bar
    
    def calc_steady_level_CL(self,W,V_inf,Sw):
        rho = self.get_density(5000)
        CL = W/(0.5*rho*V_inf*V_inf*Sw)
        return CL

    def calc_CD_from_coeff(self,CD0,CD0_L,CL, e):
        rho = self.get_density(5000)
        CD = CD0 + CD0_L*CL + ((CL*CL)/(np.pi*e*self.AR))
        return CD
    
    def get_T_req(self,CL,CD,W):
        return (CD/CL)*W

    def get_Pwr_req(self,Tr,W,CL):
        rho = self.get_density(5000)
        Pwr_req = Tr*np.sqrt((2*(W/self.Sw)/(rho*CL)))
        return Pwr_req



test_airplane = Airplane()
print(test_airplane.get_density(alt=5000))
