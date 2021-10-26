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
    return (CD/CL)*W

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
    return (Tr*np.sqrt((2*(W/Sw))/(rho*CL)))

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
    """Computes the power required for steady level flight assuming small thrust angles. The output
     units for power depend on what is input for weight, density, and wing area. 
    
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
    return (np.sqrt(2)*(term1 + term2 + term3)*W*np.sqrt((W/Sw)/rho))

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

    Returns:
        list: [minimum power required, CL that is required when flying at min power]      
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
    """Computes the velocity required for minimum power consumption given certain 
    aircraft parameters from eq 3.7.4 of Phillips "Mechanics of Flight".

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

#minimum sink airspeed
def velocity_min_sink(e,Ra,CD0,CD0_L,W,Sw,rho):
    """Computes the airspeed to fly for minimum sink or maximum endurance given certain
     aircraft parameters from eq 3.8.1 of Phillips "Mechanics of Flight".

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density
    """
    return (velocity_min_pwr(e,Ra,CD0,CD0_L,W,Sw,rho))

#maximum endurance airspeed
def velocity_max_endur(e,Ra,CD0,CD0_L,W,Sw,rho):
    """Computes the airspeed to fly for maximum endurance or minimum sink given certain
     aircraft parameters from eq 3.8.1 of Phillips "Mechanics of Flight".

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density
    """
    return (velocity_min_pwr(e,Ra,CD0,CD0_L,W,Sw,rho))

#minimum drag velocity
def velocity_min_drag(e,Ra,CD0,W,Sw,rho):
    """Computes the airspeed to fly for minimum drag given certain aircraft parameters
     from eq 3.8.2 of Phillips "Mechanics of Flight".

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density
    """
    #equation 3.8.2 from Warren Phillips "Mechanics of Flight"
    frac1 = (np.sqrt(2)/((np.pi*e*Ra*CD0)**(1/4)))
    frac2 = np.sqrt((W/Sw)/rho)
    return (frac1*frac2)

#Best Glide Airspeed
def velocity_best_glide(e,Ra,CD0,CD0_L,W,Sw,rho,headwind=0):
    """Computes the airspeed to fly for gest glide or max range given certain aircraft 
    parameters from eq 3.8.2 of Phillips "Mechanics of Flight" Returns a list of 
    [best_glide_velocity, best_glide_ratio]. Assumes the small climb/glide angle approxmiations

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        CD0_L (float): The linear coefficient in the parabolic relation for drag coefficient as a function of the lift coefficient.
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): Air density
        headwind (float): Headwind airspeed. Positive indicates a headwind. Negative indicates a tailwind.

    Returns:
        list: [Best glide airspeed, best glide ratio] 
    """
    if headwind == 0:
        #with no head or tailwind then best glide airspeed is just the minimum drag airspeed.
        V_best_glide = velocity_min_drag(e,Ra,CD0,W,Sw,rho)
        #compute the max glide ratio for zero wind Equation 3.7.19 from Phillips "Mechanics of FLight"
        numerator = np.sqrt(np.pi*e*Ra)
        denominator = (2*np.sqrt(CD0) + CD0_L*numerator)
        RG_best = numerator/denominator

    else: #There is a head or tailwind
        #set some constants that will be used in the equation of the derivative of the glide ratio w.r.t Rv from eq 3.7.25 but with CD0_L
        #Can reference "Glide_ration_deriv_with_CD0_L.png" in the directory.
        Rv_hw = (headwind/(np.sqrt((W/Sw)/rho)))
        A = CD0/2
        B = CD0_L
        C = (2/(np.pi*e*Ra))
        D = Rv_hw
        h = lambda x: (-2*A*x**5 + 3*A*D*x**4 + B*D*x*x + 2*C*x - C*D)/((A*x**4 + B*x*x + C)**2)
        Rv_guess = (4/(np.pi*e*Ra*CD0))**(1/4) #use zero wind Velocity ratio (eq. 3.7.24) as the initial guess for the newton solver 
        Rv = newton(h,Rv_guess,maxiter=10000)
        V_best_glide = Rv*np.sqrt((W/Sw)/rho)
        
        #compute the max glide ratio with wind from eq 3.7.25 but with CD0_L
        RG_best = (1 - (Rv_hw/Rv))/(A*Rv*Rv + CD0_L + C/(Rv*Rv))

    return [V_best_glide,RG_best]

#Maximum Range Airspeed
def velocity_max_range(e,Ra,CD0,W,Sw,rho):
    """Computes the airspeed to fly for max range or gest glide given certain aircraft parameters from eq 3.8.2 of Phillips "Mechanics of Flight".

    Args:
        e (float): Oswald efficiency factor (0-1)
        Ra (float): Aspect ratio of the main wing
        CD0 (float): Drag coefficient at zero lift
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): air density
    """
    return (velocity_min_drag(e,Ra,CD0,W,Sw,rho))

#Stall Speed
def velocity_stall(CL_max,W,Sw,rho):
    """Computes the stall speed of the aircraft from eq 3.8.3 of Phillips "Mechanics of Flight"

    Args:
        CL_max (float): Maximum possible lift coefficient of the aircraft
        W (float): Weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): Air density
    """
    #Equation 3.8.3 from Warren Phillips "Mechanics of Flight"
    return (np.sqrt(2/CL_max)*np.sqrt((W/Sw)/rho))


#----------------------Turns and Loads-------------------------

#Turning radius 
def turn_radius(V,phi,g,gamma=0,small_climb_angle=True):
    """Computes the turning radius of an aircraft from eq. 3.9.6 of Phillips "Mechanics of Flight"

    Args:
        V (float): airspeed
        phi (float): bank angle in radians
        g (float): acceleration due to gravity
        gamma (float): climb angle in radians (measured from the local horizontal to the flight path). Defaults to 0. Only used if small_climb_angle is set to "False"
        small_climb_angle (bool): Defaults to "True". When "True", assumes small climb angle approximation. When "False", accounts for climb angles.
    """
    #Equation 3.9.6 from Warren Phillips "Mechanics of Flight"
    #Phi cannot be 90 or -90 degrees
    if phi > -1.571 and phi < 1.5706: #These are the radian values
        #check if using small climb angle approximation
        if small_climb_angle:
            R = (V*V)/(g*np.tan(phi))
        #Use gamma if not using small climb angle approximation
        else:
            R = ((V*V*np.cos(gamma))/(g*np.tan(phi)))
    else:
        raise ValueError("Turning radius equation is undefined for a value of -90 or 90 degrees. \
            \nAlso does not accept values of phi greater than 90 and less than -90 degrees")

    return R

##Turning Climb velocity (rate of climb in a turn) computed using Thrust available and other aircraft properties
def turning_climb_rate_arspd(Ta,W,V,CD0,CD0_L,rho,Sw,e,Ra,phi,gamma=0,small_climb_angle=True):
    """Computes the turning climb rate based on the thrust available assuming small thrust angles.
    
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
        phi (float): bank angle in radians
        gamma (float): climb angle in radians (measured form the local horizontal to the flight path). Defaults to 0. Only used if small_climb_angle is set to "False"
        small_climb_angle (bool): Defaults to "True". When "True", assumes small climb angle approximation. When "False", accounts for climb angles.
    """
    if phi > -1.571 and phi < 1.5706: #These are the radian values

        if small_climb_angle: #small climb angle approximations
            #Eq. 3.9.11 from Warren Phillips "Mechanics of Flight"
            term1 = Ta*V
            term2 = 0.5*rho*V*V*V*Sw*CD0
            term3 = (W/(np.cos(phi)))
            term4 = CD0_L*V*term3
            term5 = (1/(0.5*np.pi*e*Ra*rho*V*Sw))*term3*term3
            num = term1 - (term2 + term4 + term5)
            
            Vc = num/W

        else:
            #Eq. 3.9.9 from Warren Phillips "Mechanics of Flight"
            term1 = Ta*V
            term2 = 0.5*rho*V*V*V*Sw*CD0
            term3 = ((W*np.cos(gamma))/(np.cos(phi)))
            term4 = CD0_L*V*term3
            term5 = (1/(0.5*np.pi*e*Ra*rho*V*Sw))*term3*term3
            num = term1 - (term2 + term4 + term5)
            
            Vc = num/W
    else:
        raise ValueError("Turning radius equation is undefined for a value of -90 or 90 degrees. \
            \nAlso does not accept values of phi greater than 90 and less than -90 degrees")
    
    return Vc

#Stall Limited Bank Angle
def stall_limited_bank(W,rho,V,Sw, CL_max):
    """Computes the stall limited bank angle of an aircraft (in radians) using eq. 3.9.12 from Warren Phillips "Mechanics of Flight"

    Args:
        W (float): Weight of the aircraft
        rho (float): Air density
        V (float): Airspeed
        Sw (float): Area of the main wing
        CL_max (float): Max lift coefficient of the aircraft
    """
    #Eq. 3.9.12 from Warren Phillips Mechanics of Flight
    stall_limited_phi = np.arccos(W/(0.5*rho*V*V*Sw*CL_max))
    
    return stall_limited_phi

#Load Factor Definition
def load_factor(rho,V,Sw,CL,W):
    """Computes the load factor of the aircraft.

    Args:
        rho (float): Air density
        V (float): Airspeed
        Sw (float): Area of the main wing
        CL (float): Lift coefficient during the maneuver
        W (float): Weight of the aircraft
    """
    #Eq. 3.9.13 from Warren Phillips Mechanics of Flight
    num = 0.5*rho*V*V*Sw*CL
    n = num/W #load factor

    return n
