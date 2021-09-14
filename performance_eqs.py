import numpy as np

#----------------------Thrust-------------------------

#L/D max
def L2D_max(e,Ra,CD0,CD0_L):
    """ Computes the maximum lift-to-drag ratio given an oswald efficiency, aspect ratio, and the coefficiencts of drag with respect to lift. """

    #Equation 3.2.13 from Warren Phillips "Mechanics of Flight"
    num = np.sqrt(np.pi*e*Ra)
    denom = 2*np.sqrt(CD0) + CD0_L*np.sqrt(np.pi*e*Ra)
    return (num/denom)

#Thrust Required as a function of lift and drag coefficients and weight
def thrust_req_coeff(CL,CD,W):
    """ Computes the thrust required for steady level flight assuming small thrust angles given the lift and drag coefficients and weight. 
    The units for thrust are the same as the input units used for the weight."""

    #Equation 3.2.23 from Warren Phillips "Mechanics of Flight"
    return(CD/CL)*W

#Thrust Required as a function of velocity and other aircraft parameters
def thrust_req_vel(rho,V,CD0,CD0_L,W,Sw,e,Ra):
    """ Computes the thrust required for steady level flight assumeing small thrust angles given: air density, velocity, drag coefficients, weight, wing area, oswald efficiency, 
    and aspect ratio. The units for thrust are the same as the input units used for the weight. """

    #Equation 3.2.25 from Warren Phillips "Mechanics of Flight"
    term1 = ((0.5*rho*V*V*CD0)/(W/Sw))
    term3 = ((W/Sw)/(0.5*np.pi*e*Ra*rho*V*V))
    return ((term1 + CD0_L + term3)*W)

#Minimum thrust required for small thrust angles
def thrust_req_min(CD0,CD0_L,e,Ra,W):
    """ Computes the minimum thrust required for steady level flight assuming small thrust angles given: drag coefficients, oswald efficiency, aspect ratio, and weight
    The units for thrust are the same as the input units used for the weight. """

    #Equation 3.2.26 from Warren Phillips "Mechanics of Flight"
    return ((2*np.sqrt((CD0/(np.pi*e*Ra))) + CD0_L)*W)

#----------------------Power-------------------------

#Power required as a function of thrust and aircraft params
def power_req_basic(thrust_req,W,Sw,rho,CL):
    """ Computes the power required for steady level flight assuming small thrust angles given: thrust required, weight, wing area, air density, and lift coefficient.
    The output units depend on what is input for weight, wing area, air density, and thrust required. """

    #Equation 3.3.4 from Warren Phillips "Mechanics of Flight"
    return(thrust_req*np.sqrt((2*(W/Sw))/(rho*CL)))

#Power required as a function of airspeed and aircraft params
def power_req_airspeed(CD0,CD0_L,rho,V,W,Sw,e,Ra):
    """ Computes the power required for steady level flight assuming small thrust angles given: drag coefficients, air density, velocity, weight, wing area, oswald efficiency,
    and aspect ratio. The output units for power depend on what is input for the weight, velocity, density, and wing area. """

    #Equation 3.3.7 from Warren Phillips "Mechanics of Flight"
    term1 = (CD0*rho*V**3)/(2*(W/Sw))
    term3 = (2*(W/Sw))/(np.pi*e*Ra*rho*V)
    return ((term1 + CD0_L*V + term3)*W)

####This function needs tested
#Power required as a function of CL and aircraft params
def power_req_CL(CL,CD0,CD0_L,e,Ra,W,Sw,rho):
    """ Computes the power required for steady level flight assuming small thrust angles given: Lift coefficient, drag coefficients, oswald efficiency 
    aspect ratio, weight, wing area, and air density. The output units for power depend on what is input for weight, density, and wing area. """

    #Equation 3.3.6 from Warren Phillips "Mechanics of Flight"
    term1 = (CD0/(CL**(3/2)))
    term2 = (CD0_L/np.sqrt(CL))
    term3 = (np.sqrt(CL)/(np.pi*e*Ra))
    return(np.sqrt(2)*(term1 + term2 + term3)*W*np.sqrt((W/Sw)/rho))

####This function needs tested
#Minimum Power Required for small angles
def power_req_min(CD0,CD0_L,e,Ra,W,Sw,rho):
    """ Computes the minimum power required for steady level flight assuming small thrust angles given: drag coefficients, oswald efficiency, aspect ratio 
    weight, wing area, and density. The function returns the minimum power and CL for minimum power. The output units for power depend on what is 
    input for weight, density, and wing area. """

    #Compute the CL that is required when the airplane is flying at minimum power
    #Equation 3.3.10 from Warren Phillips "Mechanics of Flight"
    term2 = np.sqrt(CD0_L**2 + ((12*CD0)/(np.pi*e*Ra)))
    CL_min_power = ((np.pi*e*Ra)/2)*(CD0_L + term2)

    #Compute the minimum power required by passing the necessary CL for minimum power to the power req_CL function
    Pr_min = power_req_CL(CL_min_power,CD0,CD0_L,e,Ra,W,Sw,rho)
    return [Pr_min,CL_min_power]
