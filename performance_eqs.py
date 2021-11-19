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

#Maneuvering Speed
def velocity_maneuver(CL_max,W_max,Sw,rho,n_pll):
    """Computes the maneuvering speed of the aircraft (the airspeed at which the airplane can perform the tightest and fastest
    maneuvers without stalling the wing or risking structural damage). Interestingly enough, maneuvering speed for an aircraft
    is also the speed for the minimum turn radius and the speed for the maximum rate of turn.

    Args:
        CL_max (float): Maximum possible lift coefficient of the aircraft
        W_max (float): Max gross weight of the aircraft
        Sw (float): Area of the main wing
        rho (float): Air density
        n_pll (float): positive load limit (defined by the maximum structural lift force/maximum gross weight.)

    Returns:
        [type]: [description]
    """
    #Equation 3.9.25 from Warren Phillips "Mechanics of Flight"
    #Compute the stall speed at maximum gross weight first
    V_smgw = velocity_stall(CL_max,W_max,Sw,rho)
    V_maneuver = np.sqrt(n_pll)*V_smgw

    return V_maneuver

#----------------------Turns and Loads-------------------------

#Turning radius from the bank angle and airspeed.
def turn_radius_from_bank_angle(V,phi,g,gamma=0,small_climb_angle=True):
    """Computes the turning radius of an aircraft in a level coordinated turn. From eq. 3.9.6 of Phillips "Mechanics of Flight"

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

#Turning radius from the load factor and airspeed.
def turn_radius_from_load_factor(V,load,g,gamma=0,small_climb_angle=True):
    """Computes the turning radius of an aircraft in a level coordinated turn. From eq. 3.9.19 and 3.9.20 of Phillips "Mechanics of Flight"
    Not valid for climbing or descending turns.

    Args:
        V (float): airspeed
        load (float): load factor (g's, defined by: Lift/Weight)
        g (float): acceleration due to gravity
        gamma (float): climb angle in radians (measured from the local horizontal to the flight path). Defaults to 0. Only used if small_climb_angle is set to "False"
        small_climb_angle (bool): Defaults to "True". When "True", assumes small climb angle approximation. When "False", accounts for climb angles.
    """
    #Equation 3.9.19 and 3.9.20 from Warren Phillips "Mechanics of Flight"
    #check if using small climb angle approximation
    if small_climb_angle:
        if load == 1: #Check to make sure you are not dividing by zero
            raise ValueError("Level turning radius equation is undefined for a value of 1 for the load factor.")
        
        else:
            R = (V*V)/(g*np.sqrt(load*load - 1))
   
    else: #Use gamma if not using small climb angle approximation
        if (load == 1 and gamma == 0): #check to make sure you are not dividing by zero
            raise ValueError("Level turning radius equation is undefined for a value of 1 for the load factor and zero climb angle.")

        else:
            R = ((V*V*(np.cos(gamma)**2))/(g*np.sqrt(load*load - (np.cos(gamma)**2))))

    return R

##Turning Climb velocity (rate of climb in a turn) computed using Thrust available and other aircraft properties
def turning_climb_rate_arspd(Ta,W,V,CD0,CD0_L,rho,Sw,e,Ra,phi,gamma=0,small_climb_angle=True):
    """Computes the turning climb rate (rate of climb in a turn) based on the thrust available assuming small thrust angles.
    
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

#Load Factor for a nearly level turn
def load_factor_level_turn(phi):
    """Approximates the load factor of an aircraft in a nearly level steady coordinated turn. usses the small climb angle approximation.

    Args:
        phi (float): Bank angle in radians
    """
    #Eq. 3.9.14 from Warren Phillips Mechanics of Flight
    n = 1/(np.cos(phi))
    
    return n

#Load limited bank angle
def load_limited_bank(W,W_max,n_pll):
    """Computes the load limited bank angle for a coordinated turn (ie. maximum bank anle that can be used without exceeding structural limits). 
    Assumes small climb angle approximation.

    Args:
        W (float): Weight of the Aircraft
        W_max (float): Max gross weight of the aircraft
        n_pll (float): positive load limit (defined by the maximum structural lift force/maximum gross weight.)

    Returns:
        float: load limited bank angle in radians
    """
    #Eq. 3.9.17 from Warren Phillips Mechanics of Flight
    load_limited_phi = np.arccos((W/(n_pll*W_max)))
    
    return load_limited_phi
 
#Stall-Limited Load Factor maximum
def max_load_factor_stall_limited(rho,V,Sw,CL_max,W):
    """Computes the maximum possible load factor that is limited by stall. Generally relevant at lower airspeeds.

    Args:
        rho (float): Air density
        V (float): Airspeed
        Sw (float): Area of the main wing
        CL_max (float): Maximum lift coefficient of the aircraft
        W (float): Weight of the aircraft
    """
    #Eq. 3.9.21 from Warren Phillips Mechanics of Flight
    n_max = load_factor(rho,V,Sw,CL_max,W)

    return n_max

#Stall-Limited Load Factor minimum
def min_load_factor_stall_limited(rho,V,Sw,CL_min,W):
    """Computes the minimum possible load factor that is limited by stall. Generally relevant at lower airspeeds.

    Args:
        rho (float): Air density
        V (float): Airspeed
        Sw (float): Area of the main wing
        CL_min (float): Minimum lif coefficient of the aircraft
        W (float): Weight of the aircraft
    """
    #Eq. 3.9.22 from Warren Phillips Mechanics of Flight
    n_min = load_factor(rho,V,Sw,CL_min,W)

    return n_min

#Structural Load Factor Maximum
def max_load_factor_structure_limited(W,W_max,n_pll):
    """Computes the maximum load factor that is limited by the structural constraints of the aircraft. 

    Args:
        W (float): Weight of the Aircraft
        W_max (float): Max gross weight of the aircraft
        n_pll (float): positive load limit (defined by the maximum structural lift force/maximum gross weight.)
    """
    #Eq. 3.9.23 from Warren Phillips Mechanics of Flight
    n_max_structure = n_pll*(W_max/W)
    
    return n_max_structure

#Structural Load Factor Minimum
def min_load_factor_structure_limited(W,W_max,n_nll):
    """Computes the minimum load factor that is limited by the structural constraints of the aircraft. 

    Args:
        W (float): Weight of the Aircraft
        W_max (float): Max gross weight of the aircraft
        n_nll (float): negative load limit (defined by the maximum negative structural lift force/maximum gross weight.)
    """
    #Eq. 3.9.23 from Warren Phillips Mechanics of Flight
    n_max_structure = n_nll*(W_max/W)
    
    return n_max_structure

#----------------------Takeoff and Landing Performance-------------------------

#Lift Off Speed
def velocity_lift_off(CL_max_config,W,Sw,rho):
    """Computes the liftoff speed of the aircraft. Defined as 1.1 times the stall speed of an airplane

    Args:
        CL_max_config (float): Maximum lift coefficient of the aircraft for take-off configuration
        W (float): Weight of the Aircraft
        Sw (float): Area of the main wing
        rho (float): Air density
    """

    #Eq. 3.10.20 from Warren Phillips Mechanics of Flight
    V_LO = 1.1*velocity_stall(CL_max_config,W,Sw,rho)

    return V_LO

#You need to program Eq. 3.10.12 but I think this will require a re-work of the equations
## as it stands Eq. 3.10.12 does not account for the air density effects on thrust, just lift and drag.
## I think I can rework the Eqs. 3.10.7 - 3.10.11 to account for the density in 3.10.7. Then code up an 
## algorithm that would discretize the starting velocity and lift off velocity into multiple smaller dV's,
## perform the "integration" and compare his results with Ex. 3.10.2

#Discrete distance calculator
def discrete_accel_decel_distance_calculator(tau, T0, T1, T2, V_start, V_end, V_hw, W, u_r, Sw, CL, CD, rho,\
    rho_0=0.0023769, a=1, g=32.17405, T_ref_vel=[]):
    """Computes the acceleration or deceleration distance of an aircraft given a starting and ending velocity. The function
    follows a modified version of Eqs. 3.10.8 - 3.10.19 in Warren Phillips "Mechanics of Flight". The modifications can 
    be found in the document listed in "Useful_reference_images". The doc is called "3.10 equations reworks.pdf".

    Args:
        tau (float): Throttle percentage. Range: 0 - 1. Where zero is 0 % throttle and one is 100% throttle
        T0 (float): First experimentally determined coefficient of thrust in the re-worked eq. (3.10.7).
        Can be a single value, or a list of values for a given range of velocities.
        T1 (float): Second experimentally determined coefficient of thrust in the re-worked eq. (3.10.7).
        Can be a single value, or a list of values for a given range of velocities.
        T2 (float): Third experimentally determined coefficient of thrust in the re-worked eq. (3.10.7).
        Can be a single value, or a list of values for a given range of velocities.
        V_start (float]): Starting velocity for distance calculation
        V_end (float): Final velocity for distance calculation
        V_hw (float): Velocity of the direct headwind component of wind
        W (float): Weight of the aircraft
        u_r (float): Coefficient of friction for rolling tires
        Sw (float): Area of the main wing
        CL (float): Lift coefficient for the configuration of the aircraft during acceleration or deceleration
        CD (float): Drag coefficient for the configuration of the aircraft during acceleration or deceleration
        rho (float): Air density
        rho_0 (float, optional): Reference air density. Defaults to 0.0023769 for English units.
        a (float): Experimentally determined exponent for density ratio (rho/rho_0). Assume 1 if unknown.
        g (float, optional): Acceleration due to gravity. Defaults to 32.17405 for English units
        T_ref_vel (list of float, optional): Thrust reference velocities to be used for interpolation if a list of values
        is given for T0, T1, and T2. All four lists must be the same dimension. Defaults to an empty list.

    Returns:
        float: The acceleration or deceleration distance for a given start and stop velocity.
    """
    #Notes: You still need to discretize the start and end velocity and create a summing loop aor keep track of the i'th 
    ##iteration. 
    ##You also need to think about relative velocities if there is a headwind. What does that mean for your start and en
    ##velocities?
    ##Good News! Your function as it stands matches what Phillips has in his Ex. 3.10.2

    #Fist thing you need to do is check if T0, T1, T2 are lists.
    ##If the are then you need to check if T0, T1, T2, and T_ref_vel are all the same length for interpolation
    placeholder = T_ref_vel
    
    #If T0, T1, T2 are not lists
    #Following Equations for 3.10.12 from Warren Phillips "Mechanics of Flight"
    #Compute K0-K2
    K0 = ((tau*T0*(rho/rho_0)**a)/W) - u_r
    K1 = (tau*T1*(rho/rho_0)**a)/W
    K2 = ((tau*T2*(rho/rho_0)**a)/W) + ((rho*Sw)/(2*W))*(CL*u_r - CD)

    #Compute KR and f values
    KR = 4*K0*K2 - K1*K1
    sr_nKR = np.sqrt(-KR)
    f1 = K0 + K1*V_start + K2*V_start*V_start
    f2 = K0 + K1*V_end + K2*V_end*V_end
    f1_p = K1 + 2*K2*V_start
    f2_p = K1 + 2*K2*V_end

    #Compute the KW and KT values
    num = (f2_p - sr_nKR)*(f1_p + sr_nKR)
    den = (f2_p + sr_nKR)*(f1_p - sr_nKR)
    KW = (1/sr_nKR)*np.log(num/den)
    KT = (1/(2*K2))*np.log(f2/f1) - ((K1*KW)/(2*K2))

    #Compute the distance
    delta_s = ((KT - V_hw*KW)/g)

    s = delta_s
    return s
#You would also need to program 3.10.25. You may also need to figure out how to rework Eqs. 3.10.21 - 3.10.24
## so they can be used in the algorithms
#Finally, I think you can also program Eq. 3.10.39 that is a less computationally expensive first approx with no wind
