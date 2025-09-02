import numpy as np

# pip install seabird (read CTD data from UNOLS)
from seabird.cnv import fCNV

# Read raw CTD profile to get sound speed values
def CTDSoundSpeed(infile):
    """
    Reads a CTD profile from a file and calculates the sound speed using the Mackenzie (1981) equation.
    """
    profile = fCNV(infile)  
    #print(profile.keys())
    latCTD=profile.attributes['LATITUDE']
    lonCTD=profile.attributes['LONGITUDE']
    datetime=profile.attributes['datetime']
    print("Profile %s \n   taken at %s, %s\n   Date Time %s" % (infile, latCTD, lonCTD, datetime))
    # from Mackenzie (1981)
    C_Mac = 1448.6 + 4.618*profile['TEMP2'] + \
        - 0.0541*profile['TEMP2']**2 + \
        + 0.0002*profile['TEMP2']**3 + \
        + 1.34*(profile['PSAL2'] - 35) + \
        + 0.016*profile['DEPTH']
    return  profile, C_Mac
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
    
