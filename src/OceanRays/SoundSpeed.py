import numpy as np


# Read raw CTD profile to get sound speed values
def CTDSoundSpeed(profile,model='Mackenzie'):
    """
    Reads a CTD profile from a file and calculates the sound speed using the Mackenzie (1981) equation.
    """
    def C_Mackenzie(T, S, D):
        """
        Calculate sound speed using Mackenzie (1981) equation.
        
        Args:
            T (float or array): Temperature in degrees Celsius.
            S (float or array): Salinity in PSU.
            D (float or array): Depth in meters.
        
        Returns:
            float or array: Sound speed in m/s.
        """
        C = 1448.96 + 4.591*T - 0.05304*T**2 + 0.0002374*T**3 + 1.34*(S - 35) + 0.0163*D +1.675e-7*D**2 - 0.01025*T*(S - 35) - 7.139e-13*T*D**3
        return C
    def C_Leroy(T, S, D,L=45):
        """
        Calculate sound speed using Leroy (2008) equation.
        implementation from: https://gorbatschow.github.io/SonarDocs/sound_speed_sea_npl.en/
        Args:
            T (float or array): Temperature in degrees Celsius.
            S (float or array): Salinity in PSU.
            D (float or array): Depth in meters.
            L (float): Latitude in degrees. Default is 45.
        
        Returns:
            float or array: Sound speed in m/s.
        """
        C = 1402.5 + 5*T - (5.44e-2)*(T**2) + (2.1e-4)*(T**3) + \
            1.33*S - (1.23e-2)*S*T + (8.7e-5)*S*(T**2) + \
            (1.56e-2)*D + (2.55e-7)*(D**2) - (7.3e-12)*(D**3) + \
            (1.2e-6)*D*(L-45) - (9.5e-13)*T*(D**3) + (3e-7)*(T**2)*D + \
            (1.43e-5)*S*D
        return C
    def C_Coppens(T, S, D):
        """
        Calculate sound speed using Coppens (1981) equation.
        implentation from: https://gorbatschow.github.io/SonarDocs/sound_speed_sea_coppens.en/
        
        Args:
            T (float or array): Temperature in degrees Celsius.
            S (float or array): Salinity in PSU.
            D (float or array): Depth in meters.
        
        Returns:
            float or array: Sound speed in m/s.
        """
        d = D*(1e-3);
        t = T*(1e-1);

        C = 1449.05 + 45.7*t - 5.21*(t**2) + 0.23*(t**3) + \
                + (1.333 - 0.126*t + 0.009*(t**2))*(S - 35) + \
                + (16.23 + 0.253*t)*d + (0.213-0.1*t)*(d**2) + \
                + (0.016 + 0.0002*(S-35))*(S-35)*t*d;

        return C
# Below are considerably different.  Pressure may be wrong scale.
    def C_DelGrosso(T, S, D, P=None):
        """
        Calculate sound speed using Del Grosso (1974) equation.
        
        Args:
            T (float or array): Temperature in degrees Celsius.
            S (float or array): Salinity in PSU.
            D (float or array): Depth in meters.
        Using implementation from: https://gorbatschow.github.io/SonarDocs/sound_speed_sea_delgrosso.en/
        Returns:
            float or array: Sound speed in m/s.
        """
        C000 = 1402.392
        CT1 = 5.012285
        CT2 = -0.0551184
        CT3 = 0.221649e-3
        CS1 = 1.329530
        CS2 = 0.1288598e-3
        CP1 = 0.1560592
        CP2 = 0.2449993e-4
        CP3 = -0.8833959e-8
        CST = -0.01275936
        CTP = 0.6353509e-2
        CT2P2 = 0.2656174e-7
        CTP2 = -0.1593895e-5
        CTP3 = 0.5222483e-9
        CT3P = -0.4383615e-6
        CS2P2 = -0.1616745e-8
        CST2 = 0.9688441e-4
        CS2TP = 0.4857614e-5
        CSTP =-0.3406824e-3
        
        #if P is None:
        #    P=D/10; # pressure in MPa (1 m = 0.1 MPa)
        #p = P*0.0980665; # from dbar to atm

        if P is None:
            p = D*0.1 # from dbar to atm
        else:
            p = P*.1 

        CT = CT1*T + CT2*(T**2) + CT3*(T**3)
        CS = CS1*S + CS2*(S**2)
        CP = CP1*p + CP2*(p**2) + CP3*(p**3)
        CSTP = CTP*T*p + CT3P*(T**3)*p + CTP2*T*(p**2) + \
            + CT2P2*(T**2)*(p**2) + CTP3*T*(p**3) + \
            + CST*S*T + CST2*S*(T**2) + CSTP*S*T*p + \
            + CS2TP*(S**2)*T*p + CS2P2*(S**2)*(p**2)
        C = C000 + CT + CS + CP + CSTP
        return C    
    def C_Chen_Millero(T, S, D, P=None):
        """
        Calculate sound speed using Chen and Millero (1977) equation.
        
        Args:
            T (float or array): Temperature in degrees Celsius.
            S (float or array): Salinity in PSU.
            D (float or array): Depth in meters.
        Implementation from : https://gorbatschow.github.io/SonarDocs/sound_speed_sea_unesco.en/
        Returns:
            float or array: Sound speed in m/s.
        """
        C00=1402.388; C01=5.03830; C02=-5.81090e-2; C03=3.3432e-4 
        C04=-1.47797e-6; C05=3.1419e-9
        C10=0.153563; C11=6.8999e-4; C12=-8.1829e-6; C13=1.3632e-7
        C14=-6.1260e-10; C20=3.1260e-5; C21=-1.7111e-6
        C22=2.5986e-8; C23=-2.5353e-10; C24=1.0415e-12; C30=-9.7729e-9
        C31=3.8513e-10; C32=-2.3654e-12
        A00=1.389; A01=-1.262e-2; A02=7.166e-5; A03=2.008e-6
        A04=-3.21e-8
        A10=9.4742e-5; A11=-1.2583e-5; A12=-6.4928e-8; A13=1.0515e-8
        A14=-2.0142e-10
        A20=-3.9064e-7; A21=9.1061e-9; A22=-1.6009e-10; A23=7.994e-12
        A30=1.100e-10; A31=6.651e-12; A32=-3.391e-13
        B00=-1.922e-2; B01=-4.42e-5; B10=7.3637e-5; B11=1.7950e-7
        D00=1.727e-3
        D10=-7.9836e-6

        if P is None:
            p = D*0.1  # from dbar to atm
        else:
            p = P*.1
        
        Cw1 = (C00 + C01*T + C02*(T**2) + C03*(T**3) + C04*(T**4) + C05*(T**5)
              + (C10 + C11*T + C12*(T**2) + C13*(T**3) + C14*(T**4))*p
              + (C20 + C21*T + C22*(T**2) + C23*(T**3) + C24*(T**4))*(p**2)
              + (C30 + C31*T + C32*(T**2))*(p**3))

        A1 = (A00 + A01*T + A02*(T**2) + A03*(T**3) + A04*(T**4)
             + (A10 + A11*T + A12*(T**2) + A13*(T**3) + A14*(T**4))*p
             + (A20 + A21*T + A22*(T**2) + A23*(T**3))*(p**2)
             + (A30 + A31*T + A32*(T**2))*(p**3))

        B1 = B00 + B01*T + (B10 + B11*T)*p

        D1 = D00 + (D10*p)
        
        C = Cw1 + A1*S + B1*(S**(3/2)) + D1*(S**2)
        return C

    #profile = fCNV(infile)  
    #print(profile.keys())
    latCTD=profile.attributes['LATITUDE']
    lonCTD=profile.attributes['LONGITUDE']
    datetime=profile.attributes['datetime']
    #print("Profile:\n taken at %s, %s\n   Date Time %s" % (latCTD, lonCTD, datetime))
    # from Mackenzie (1981)
    if model=='Mackenzie':
        C = C_Mackenzie(profile['TEMP'], profile['PSAL'], profile['DEPTH'])
    elif model=='Leroy':
        C = C_Leroy(profile['TEMP'], profile['PSAL'], profile['DEPTH'],latCTD)
    elif model=='Coppens':
        C = C_Coppens(profile['TEMP'], profile['PSAL'], profile['DEPTH'])   
    elif model=='DelGrosso':
        C = C_DelGrosso(profile['TEMP'], profile['PSAL'], profile['DEPTH'], profile['PRES']) 
    elif model=='Chen_Millero':
        C = C_Chen_Millero(profile['TEMP'], profile['PSAL'], profile['DEPTH'], profile['PRES'])
    else:
        raise ValueError("Model not recognized. Use one of the recognized models ('Mackenzie', 'Leroy', 'Coppens', 'DelGrosso', or 'Chen_Millero').")
    return C
# Sound Speed Gradient Calculation with interpolation along z
def SoundSpeedGrad(z, zProfile, cProfile, ssp_only=False):
    """
    Calculates the interpolated sound speed and its gradient at a specific depth 
    using linear interpolation.

    Args:
        z (float): Depth at which to calculate sound speed (m).
        zProfile (array): Array of depths from the CTD profile (m).
        cProfile (array): Array of sound speeds corresponding to the depths (m/s).
        ssp_only (bool): If True, return only the sound speed profile.

    Returns:
        tuple: Interpolated sound speed and its gradient (if not ssp_only).
    """
    # Interpolate to desired depths
    zMin=zProfile.min()
    zMax=zProfile.max()
    if z < zMin:
        #print("z is below zMin.  Extrapolating sound speed profile using gradient of nearest reported values")
        cInterpP1 = np.interp(zMin+1, zProfile, cProfile) # plus 1 m
        cInterp0 = np.interp(zMin, zProfile, cProfile) # at 1st point
        zdiff=zMin-z
        cInterpGrad = (cInterpP1 - cInterp0)
        cInterp=cInterp0 + zdiff*cInterpGrad
    elif z > zMax:
        #print("z is above zMax.  Extrapolating sound speed profile using gradient of nearest reported values")
        cInterpM1 = np.interp(zMax-1, zProfile, cProfile) # minus 1 m
        cInterpLast = np.interp(zMax, zProfile, cProfile) # last point
        zdiff=z-zMax
        cInterpGrad = (cInterpLast - cInterpM1)
        cInterp=cInterpLast + zdiff*cInterpGrad
    else:
        cInterpM1 = np.interp(z-1, zProfile, cProfile) # minus 1 m
        cInterpP1 = np.interp(z+1, zProfile, cProfile) # plus 1 m
        cInterpGrad = (cInterpP1 - cInterpM1) / 2.0         
        cInterp = np.interp(z, zProfile, cProfile) 
    
    if ssp_only:
        return cInterp
    else:
        return cInterp, cInterpGrad
# Theoretical Sound Speed Profile and Gradient
def MunkSoundSpeedGrad(z,c0=1500, z_min=1300, B=1300, epsilon=0.0057, ssp_only=False):
    """
    Defines the Munk sound speed profile and gradiant.
    The gradient is essential for determining the curvature of the sound rays.

    Args:
        z (float): Depth in meters (positive downwards).

        # Munk profile parameters
        c0 = Reference sound speed (m/s)
        z_min =  Depth of the sound channel minimum (m)

        B = Scale depth, related to the channel's width (m)
        epsilon =  Dimensionless parameter

        ssp_only (bool): If True, return only the sound speed profile.

    Returns:
        float: Sound speed in m/s.
        float: The gradient of sound speed with respect to depth (1/s).
    """
    eta = 2 * (z - z_min) / B
    grad = c0 * epsilon * (2 / B) * (1 - np.exp(-eta))
    c_profile = c0 * (1 + epsilon * (eta + np.exp(-eta) - 1))
    if ssp_only :
      return c_profile
    else:
      return c_profile, grad
    
