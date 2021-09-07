import numpy as np
import std_atmos_Eng
import std_atmos_SI

class Airplane():
    def __init__(self):
        self.wingspan = 10  #ft
        self.c_bar = 2 #ft
        self.T0 = 5000 #lbf
        self.AR = self.wingspan/self.c_bar
        self.e = 0.8 #Oswald efficiency
        self.Sw = self.wingspan*self.c_bar #wing area
        self.units = 'English'

    def get_density(self,alt):

        #check which units you are using
        if self.units == 'English':
            rho = std_atmos_Eng.interp(alt)[3] #return the density from the interpolation 
        else: #Use SI
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
print(test_airplane.get_AR())

print(test_airplane.calc_CD_from_coeff(0.03,0.001,0.7,0.8))
print('The thrust required is: ',test_airplane.get_T_req(0.7,0.06969,2000))
print('The power required is: ',test_airplane.get_Pwr_req(199.114,2000,0.7))