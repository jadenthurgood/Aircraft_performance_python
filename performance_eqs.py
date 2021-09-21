import numpy as np
from scipy.optimize import newton

#----------------------Thrust-------------------------

#L/D max
def L2D_max(e,Ra,CD0,CD0_L):
    """ Computes the maximum lift-to-drag ratio.
    
    Args:
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing 
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.        
    """

    #Equation 3.2.13 from Warren Phillips "Mechanics of Flight"
    num = np.sqrt(np.pi*e*Ra)
    denom = 2*np.sqrt(CD0) + CD0_L*np.sqrt(np.pi*e*Ra)
    return (num/denom)

#Thrust Required as a function of lift and drag coefficients and weight
def thrust_req_coeff(CL,CD,W):
    """ Computes the thrust required for steady level flight assuming small thrust angles given the lift and drag coefficients and weight. 
    The units for thrust are the same as the input units used for the weight.
    
    Args:
        CL (float): Lift Coefficient at steady level flight
        CD (float): Drag Coefficient under the same conditions as the lift coefficient
        W (float): Weight of the aircraft
    """

    #Equation 3.2.23 from Warren Phillips "Mechanics of Flight"
    return(CD/CL)*W

#Thrust Required as a function of velocity and other aircraft parameters
def thrust_req_arspd(rho,V,CD0,CD0_L,W,Sw,e,Ra):
    """ Computes the thrust required for steady level flight assumeing small thrust angles. The units for thrust are 
    the same as the input units used for the weight. 
    
    Args:
        rho (float): Air density
        V (float): Aispeed
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.      
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing 
    """

    #Equation 3.2.25 from Warren Phillips "Mechanics of Flight"
    term1 = ((0.5*rho*V*V*CD0)/(W/Sw))
    term3 = ((W/Sw)/(0.5*np.pi*e*Ra*rho*V*V))
    return ((term1 + CD0_L + term3)*W)

#Minimum thrust required for small thrust angles
def thrust_req_min(CD0,CD0_L,e,Ra,W):
    """ Computes the minimum thrust required for steady level flight assuming small thrust angles. The units for thrust 
    are the same as the input units used for the weight. 
    
    Args:
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing
        W (float): Weight of the aircraft        
    """

    #Equation 3.2.26 from Warren Phillips "Mechanics of Flight"
    return ((2*np.sqrt((CD0/(np.pi*e*Ra))) + CD0_L)*W)

#----------------------Power-------------------------

#Power required as a function of thrust and aircraft params
def power_req_basic(Tr,W,Sw,rho,CL):
    """ Computes the power required for steady level flight assuming small thrust angles. The output units depend on what is 
    input for weight, wing area, air density, and thrust required. 
    
    Args:
        Tr (float): Thrust required for steady level flight. Can be computed using thrust_req_coeff() or thrust_req_arspd().
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): Air density
        CL (float): Lift coefficient for the given thrust required
    """

    #Equation 3.3.4 from Warren Phillips "Mechanics of Flight"
    return(Tr*np.sqrt((2*(W/Sw))/(rho*CL)))

#Power required as a function of airspeed and aircraft params
def power_req_arspd(CD0,CD0_L,rho,V,W,Sw,e,Ra):
    """ Computes the power required for steady level flight assuming small thrust angles. The output units for power depend on
    what is input for the weight, velocity, density, and wing area. 
    
    Args:
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        rho (float): Air density
        V (float): Aispeed
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing       
    """

    #Equation 3.3.7 from Warren Phillips "Mechanics of Flight"
    term1 = (CD0*rho*V**3)/(2*(W/Sw))
    term3 = (2*(W/Sw))/(np.pi*e*Ra*rho*V)
    return ((term1 + CD0_L*V + term3)*W)

####This function needs tested
#Power required as a function of CL and aircraft params
def power_req_CL(CL,CD0,CD0_L,e,Ra,W,Sw,rho):
    """Computes the power required for steady level flight assuming small thrust angles. The output units for power depend on what is input for weight, density, and wing area. 
    
    Args:
        CL (float): Lift coefficient
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing      
        rho (float): air density       
    """

    #Equation 3.3.6 from Warren Phillips "Mechanics of Flight"
    term1 = (CD0/(CL**(3/2)))
    term2 = (CD0_L/np.sqrt(CL))
    term3 = (np.sqrt(CL)/(np.pi*e*Ra))
    return(np.sqrt(2)*(term1 + term2 + term3)*W*np.sqrt((W/Sw)/rho))

####This function needs tested
#Minimum Power Required for small angles
def power_req_min(CD0,CD0_L,e,Ra,W,Sw,rho):
    """Computes the minimum power required for steady level flight assuming small thrust angles. The function returns the minimum power 
    and CL for minimum power as a list [min_pwr, CL_for_min_pwr]. The output units for power depend on what is input for weight, density, and wing area. 
    
    Args:
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density        
    """

    #Compute the CL that is required when the airplane is flying at minimum power
    #Equation 3.3.10 from Warren Phillips "Mechanics of Flight"
    term2 = np.sqrt(CD0_L**2 + ((12*CD0)/(np.pi*e*Ra)))
    CL_min_power = ((np.pi*e*Ra)/2)*(CD0_L + term2)

    #Compute the minimum power required by passing the necessary CL for minimum power to the power req_CL function
    Pr_min = power_req_CL(CL_min_power,CD0,CD0_L,e,Ra,W,Sw,rho)
    return [Pr_min,CL_min_power]

#----------------------Velocities-------------------------

####This function needs tested
#Climb velocity (rate of climb) computed using power available and power required
def climb_rate_power(pwr_avail,pwr_req,W):
    """Computes the climb rate based on the excess power available over the power required for steady level flight. 
    The function assumes small climb angles.
    
    Args:
        pwr_avail (float): Power available
        pwr_req (float): Power required for steady level flight
        W (float): Weight of the aircraft
    """
    
    #Equation 3.4.8 from Warren Phillips "Mechanics of Flight"
    return ((pwr_avail - pwr_req)/W)

#Climb velocity (rate of climb) computed using Thrust available and other aircraft properties
def climb_rate_arspd(Ta,W,V,CD0,CD0_L,rho,Sw,e,Ra,small_climb_angle=True):
    """Computes the climb rate based on the thrust available assuming small thrust angles and small climb angles.
    
    Args:
        Ta (float): Oswald efficiency factor (0-1)
        W (float): Weight of the aircraft
        V (float): Airspeed
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        rho (float): air density
        Sw (float): Area of the main wing
        e (float): Oswald efficiency (0-1)
        Ra (float): Aspect ratio of the main wing
        small_climb_angle (bool): Defaults to "True". When "True", assumes small climb angle approximation. When "False", accounts for climb angles, and sets CD0_L = 0 for computational purposes. Results are valid only when CD0_L is approx "0"
    """
    if small_climb_angle: #small climb angle approximations
        #Page 283 from Warren Phillips "Mechanics of Flight"
        term1 = (Ta/W)*V
        terms2_4 = (power_req_arspd(CD0,CD0_L,rho,V,W,Sw,e,Ra))/W
        Vc = (term1 - terms2_4)

    else: #without using the small climb angle approximation 
        #compute CL Equations come from the image titled "CL_equations_climb_angle" 
        #Comes from trying to solve the equation on page 283 of Phillips "Mechanics of Flight"
        A = ((0.5*rho*V*V*Sw)/W)
        B = (1/(np.pi*e*Ra))
        C = Ta/W

        const1 = (np.sqrt(4*A*A*B*CD0 + A*A - 4*A*B*C + 4*B*B)/(A*B*B))
        const2 = ((2*C)/(A*B))
        const3 = 1/(B*B)
        const4 = ((2*CD0)/B)

        CL1 = -((np.sqrt(-const1 + const2 - const3 - const4))/np.sqrt(2))
        CL2 = ((np.sqrt(-const1 + const2 - const3 - const4))/np.sqrt(2))
        CL3 = -np.sqrt(0.5*(const1 + const2 - const3 - const4))
        CL4 = np.sqrt(0.5*(const1 + const2 - const3 - const4))
        
        #find the CL that is positive and real.
        #put the computed CL values in a list
        CL_list = [CL1,CL2,CL3,CL4]
        #remove the NaN's from the list
        real_CL_list = [x for x in CL_list if np.isnan(x) == False]
        #select the positive value which is the only element in the list at this point. So index at [0]
        positive_CL = [y for y in real_CL_list if y >= 0][0]
        
        #compute gamma
        gamma = np.arccos(A*positive_CL)

        #compute Vc
        Vc = V*np.sin(gamma)

    return Vc

####This function needs tested
#Minimum Power Required airspeed
def velocity_min_pwr(e,Ra,CD0,CD0_L,W,Sw,rho):
    """Computes the velocity required for minimum power consumption given certain aircraft parameters from eq 3.7.4 of Phillips "Mechanics of Flight".

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density
    """
    #Equation 3.7.4 from Warren Phillips "Mechanics of Flight"
    term1 = np.pi*e*Ra*CD0_L
    term2 = np.sqrt((term1**2) + 12*np.pi*e*Ra*CD0)
    V_mdv = (2/(np.sqrt(term1 + term2)))*(np.sqrt((W/Sw)/rho))
    return V_mdv

####This function needs tested
#Sink rate given power required
def sink_rate_pwr(pwr_req,W):
    """Computes the sink rate of the aircraft from eq 3.7.3 of Phillips "Mechanics of Flight".

    Args:
        pwr_req (float): Power required for level flight. Can be computed using power_req_basic(), power_req_arspd(), or power_req_CL()
        W (float): Aircraft weight
    """
    #Equation 3.7.3 from Warren Phillips "Mechanics of Flight"
    return (pwr_req/W)

####This function needs tested
#Sink rate given other aircraft parameters
def sink_rate_arspd(CD0,CD0_L,rho,V,W,Sw,e,Ra):
    """Computes the sink rate of the aircraft from eq 3.7.3 of Phillips "Mechanics of Flight"

    Args:
        CD0 (float): Drag Coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        rho (float): Air density
        V (float): Airspeed
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect Ratio of the main wing
    """
    
    #Equation 3.7.3 from Warren Phillips "Mechanics of Flight"
    Vs = (power_req_arspd(CD0,CD0_L,rho,V,W,Sw,e,Ra))/W
    return Vs

####This function needs tested
#Minimum sink rate 
def sink_rate_min(CD0,rho,W,Sw,e,Ra):
    """Computes the minimum sink rate if the pilot flies at minimum power velocity.

    Args:
        CD0 (float): Drag Coefficient at zero lift
        rho (float): Air density
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect Ratio of the main wing
    """

    #Equation 3.7.5 from Warren Phillips "Mechanics of Flight"
    frac1 = (4*np.sqrt(2)*(CD0**(1/4)))/((3*np.pi*e*Ra)**(3/4))
    frac2 = np.sqrt((W/Sw)/rho)
    return (frac1*frac2)

#####Best Glide Airspeeds and Glide Ratios#####
print (climb_rate_arspd(6500,20000,528,0.023,0.0,0.0023769,320,0.82,9.1125,False))
