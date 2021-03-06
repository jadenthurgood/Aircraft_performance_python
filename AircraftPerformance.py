import numpy as np
from Standard_Atmosphere import std_atmos_Eng
from Standard_Atmosphere import std_atmos_SI 
import warnings

class Airplane():
    def __init__(self):
        #Aircraft Geometry
        self.wingspan = 25  #ft
        self.c_bar = 5 #ft
        self.Sw = self.wingspan*self.c_bar #wing area

        #Aricraft Characteristics
        self.CD0 = 0.03
        self.CD0_L = 0.001
        self.T0 = 5000 #lbf
        self.e = 0.8 #Oswald efficiency
        self.weight = 800
        
        #Conditions
        self.units = 'English'
        self.V_inf = 100
        self.altitude = 5000

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

    def get_AR(self,**kwargs):
        """ This function computes the aspect ratio of the aircraft. It can do this using the input file values, or the user may set the wingspan
        and average chord length as keyword arguments 'span' and 'c_bar'. If the user specifies the span and chord length in the function call it will overite the 
        values from the input file if one is given. """

        if kwargs: #check if the user input any kwargs for the aspect ratio calculation
            if ('span' in kwargs) and ('c_bar' not in kwargs): #check if the user only input a span value
               AR = kwargs['span']/self.c_bar 
            elif ('c_bar' in kwargs) and ('span' not in kwargs): #check if the user only input the average chord
                AR = self.wingspan/kwargs['c_bar']
            elif ('span' and 'c_bar' in kwargs): #Check if the user input both the span and average chord
                AR = kwargs['span']/kwargs['c_bar']
            else: #if the user has input kwargs that are neither 'span' or 'c_bar' raise an error. They are trying something incorrect
                raise ValueError("The only keyword arguments used in get_AR are 'span' and/or 'c_bar'.")

        else: #use the values from the input file to compute the aspect ration
            AR = self.wingspan/self.c_bar
        return AR
    
    def calc_steady_level_CL(self,**kwargs): #function uses the small thrust angle approximation from phillips "Mechanics of Flight" Eq (3.2.24)
        """ This function calculates CL for steady level flight assuming small thrust angles. It can do this using the input file values or the user can
        manually change the conditions using keyword arguments 'weight', 'V_inf' for freestream velocity, 'Sw' for wing area, 'alt' for altitude and 'units'.
        Any kewword input values will override the input file values. """

        #unpack the values from kwargs or from Class self
        if kwargs: #check if the user input any kwargs for stead CL calculation
            
            if 'weight' in kwargs: #if the user assigned a weight in the function call, then use it. 
                w = kwargs['weight']
            else: #Otherwise, default to the input file value
                w = self.weight

            if 'V_inf' in kwargs: #if the user assigned a freestream velocity in the function call, then use it. 
                V = kwargs['V_inf']
            else: #Otherwise, default to the input file value
                V = self.V_inf

            if 'Sw' in kwargs: #if the user assigned a reference area in the function call, then use it. 
                area_ref = kwargs['Sw']
            else: #Otherwise, default to the input file value
                area_ref = self.Sw
            
            if 'alt' in kwargs: #if the user assigned an altitude in the function call, then use it. 
                h = kwargs['alt']
            else: #Otherwise, default to the input file value
                h = self.altitude
            
            if 'units' in kwargs: #if the user assigned a units designation in the function call, then use it. 
                local_units = kwargs['units']
            else: #Otherwise, use the input file value
                pass
            
            #if the user inputs a kwarg that is not used without any kwargs that are used. Throw an errorgit 
            if (any(key in kwargs for key in ('weight','V_inf','Sw','alt','units'))): 
                None
            else:
                raise ValueError("The only keyword arguments used in calc_steady_level_CL are any combination of 'weight', 'V_inf', 'Sw', 'alt', and 'units'.")
        
        else: #just use the values from the input file
            w = self.weight
            V = self.V_inf
            area_ref = self.Sw
            h = self.altitude
        
        #Get the density
        if 'local_units' in locals():
            rho = self.get_density(h,units=local_units) #get the density with the user specified units 
        else: 
            rho = self.get_density(h)

        #compute CL
        CL = w/(0.5*rho*V*V*area_ref) #from Eq. (3.2.24)

        return CL

    def calc_CD_from_coeff(self,CL,**kwargs): #Equation from Warren Phillips "Mechanics of Flight" Eq. (3.2.7)
        """ This function calculates CD from the coefficients of the parabolic function with respect to a user input CL. The user may specify these values
        using the input file or at the time of the function call using the kwargs 'CD0', 'CD0_L', 'e' for oswald efficiency, 'AR' for aspect ratio.
        Any inputs during the time of the function call will overwrite the values from the input file. """
        
        #unpack kwargs
        if kwargs: #check if the user input any kwargs for stead CL calculation
            
            if 'CD0' in kwargs: #if the user assigned a CD0 in the function call, then use it. 
                CD0 = kwargs['CD0']
            else: #Otherwise, default to the input file value
                CD0 = self.CD0
            
            if 'CD0_L' in kwargs: #if the user assigned a CD0_L in the function call, then use it. 
                CD0_L = kwargs['CD0_L']
            else: #Otherwise, default to the input file value
                CD0_L = self.CD0_L
            
            if 'e' in kwargs: #if the user assigned a value for 'e' in the function call, then use it. 
                e = kwargs['e']
            else: #Otherwise, default to the input file value
                e = self.e
            
            if 'AR' in kwargs: #if the user assigned a value for 'AR' in the function call, then use it. 
                AR = kwargs['AR']
            else: #Otherwise, default to the input file value
                AR = self.get_AR()

            #if the user inputs a kwarg that is not used without any kwargs that are used. Throw an errorgit 
            if (any(key in kwargs for key in ('CD0','CD0_L','e','AR'))): 
                None
            else:
                raise ValueError("The only keyword arguments used in calc_CD_from_coeff are any combination of 'CD0', 'CD0_L', 'e', and 'AR'.")
        
        else: #use the values from the input file if no kwargs are given
            CD0 = self.CD0
            CD0_L = self.CD0_L
            e = self.e
            AR = self.get_AR()
        
        #compute CD
        CD = CD0 + CD0_L*CL + ((CL*CL)/(np.pi*e*AR)) #Variation from Eq (3.2.7) from Mechanics of Flight
        return CD
    
    def get_T_req(self,type='velocity',minimum=False,**kwargs): #compute the thrust required. From Phillips "Mechanics of Flight" Eqs.(3.2.22 - 3.2.26)
        """ Computes the thrust required for steady level flight. "type" determines if thrust required is a function of 'velocity', 'coefficients', or 'L/D'.
        The user may also specify if they want the the minimum thrust required by setting 'minimum' = True. """
        #Thrust required as a function of velocity
        if type == 'velocity':
            #unpack the values from kwargs or from Class self
            if kwargs: #check if the user input any kwargs for stead CL calculation
                
                if 'weight' in kwargs: #if the user assigned a weight in the function call, then use it. 
                    w = kwargs['weight']
                else: #Otherwise, default to the input file value
                    w = self.weight

                if 'V_inf' in kwargs: #if the user assigned a freestream velocity in the function call, then use it. 
                    V = kwargs['V_inf']
                else: #Otherwise, default to the input file value
                    V = self.V_inf

                if 'Sw' in kwargs: #if the user assigned a reference area in the function call, then use it. 
                    area_ref = kwargs['Sw']
                else: #Otherwise, default to the input file value
                    area_ref = self.Sw
                
                if 'alt' in kwargs: #if the user assigned an altitude in the function call, then use it. 
                    h = kwargs['alt']
                else: #Otherwise, default to the input file value
                    h = self.altitude
                
                if 'CD0' in kwargs: #if the user assigned a CD0 in the function call, then use it. 
                    CD0 = kwargs['CD0']
                else: #Otherwise, default to the input file value
                    CD0 = self.CD0
                
                if 'CD0_L' in kwargs: #if the user assigned a CD0_L in the function call, then use it. 
                    CD0_L = kwargs['CD0_L']
                else: #Otherwise, default to the input file value
                    CD0_L = self.CD0_L
                
                if 'e' in kwargs: #if the user assigned a value for 'e' in the function call, then use it. 
                    e = kwargs['e']
                else: #Otherwise, default to the input file value
                    e = self.e
                
                if 'AR' in kwargs: #if the user assigned a value for 'AR' in the function call, then use it. 
                    AR = kwargs['AR']
                else: #Otherwise, default to the input file value
                    AR = self.get_AR()
                    
                if 'units' in kwargs: #if the user assigned a units designation in the function call, then use it. 
                    local_units = kwargs['units']
                else: #Otherwise, use the input file value
                    pass
                
                #if the user inputs a kwarg that is not used without any kwargs that are used. Throw an errorgit 
                if (any(key in kwargs for key in ('weight','V_inf','Sw','alt','CD0','CD0_L','e','AR','units'))): 
                    None
                else:
                    raise ValueError("The only keyword arguments used in calc_steady_level_CL are any combination of 'weight','V_inf','Sw','alt','CD0','CD0_L','e','AR','units'.")
            
            else: #just use the values from the input file
                w = self.weight
                V = self.V_inf
                area_ref = self.Sw
                h = self.altitude
                CD0 = self.CD0
                CD0_L = self.CD0_L
                e = self.e
                AR = self.get_AR()

            ######What if we handled units with a set_units function. Tather than passing the units everytime we need it. The user could call set_units() and set the global units to 
            ######what he needed it to be before making a call to another function
            
            #Get the density
            if 'local_units' in locals():
                rho = self.get_density(h,units=local_units) #get the density with the user specified units 
            else: 
                rho = self.get_density(h)
        
            term1 = ((0.5*rho*V*V*CD0)/(w/area_ref))
            term3 = ((w/area_ref)/(0.5*np.pi*e*AR*rho*V*V))
            T_req = (term1 + CD0_L + term3)*w #Eq. (3.2.25) Mechanics of Flight     
        
        return T_req

    def get_Pwr_req(self,Tr,W,CL):
        rho = self.get_density(5000)
        Pwr_req = Tr*np.sqrt((2*(W/self.Sw)/(rho*CL)))
        return Pwr_req

test_airplane = Airplane()
print(test_airplane.get_T_req())
