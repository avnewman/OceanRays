import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pickle

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
# Ray Trace either a given sound speed profile or the Munk profile
def RayTrace(initial_angle_deg, max_range_km, initial_depth, dt=0.1,zProfile=None,cProfile=None):
    """
    Calculates the path of a single acoustic ray and its travel time to the surface.

    Args:
        initial_angle_deg (float): Initial angle of the ray in degrees from the horizontal.
        max_range_km (float): Maximum horizontal distance to trace the ray in km.
        initial_depth (float): Starting depth of the ray in meters.
        dt (float): Time step for the integration in seconds.

    Returns:
        tuple: A tuple containing (ranges, depths, travel_time).
               travel_time is the time in seconds to reach the surface, or None if it doesn't.
    """
    # Convert initial angle to radians
    theta = np.deg2rad(initial_angle_deg)

    # Initial conditions
    r = 0  # Range in meters
    z = initial_depth
    t = 0  # Time in seconds

    # Store path for plotting
    ranges = [r / 1000]
    depths = [z]

    max_range_m = max_range_km * 1000
    travel_time_to_surface = None

    # Main loop for ray tracing
    while r < max_range_m:
        if zProfile is not None and cProfile is not None:
            # using sound speed profile
            c, dc_dz = SoundSpeedGrad(z, zProfile, cProfile)
        else:
            # using Munk profile
            c, dc_dz = MunkSoundSpeedGrad(z)

        # Update ray parameters using the ray tracing equations
        r_step = c * np.cos(theta) * dt
        z_step = c * np.sin(theta) * dt
        
        # Check if the ray crosses the surface in the next step
        if z + z_step <= 0:
            # Interpolate to find the exact time and position of surface crossing
            frac = -z / z_step
            t += frac * dt
            r += frac * r_step
            z = 0
            travel_time_to_surface = t
            ranges.append(r)
            depths.append(z)
            break # Stop tracing once the surface is reached

        r += r_step
        z += z_step
        theta -= (dc_dz / c) * c * np.cos(theta) * dt
        t += dt

        ranges.append(r)
        depths.append(z)

        # Stop if the ray hits the surface (e.g., 0m)
        if z < 0:
            break

    return ranges, depths, travel_time_to_surface

def plotProfileRaypaths(ray_paths, zProfile=None, cProfile=None, sspDepths=None):
    """
    Plots the sound speed profile and the ray paths.

    Args:
        ray_paths (dict): A dictionary where keys are angles and values are (path, travel_time).
        sspDepths (np.array): An array of depth values to plot the sound speed profile.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8), gridspec_kw={'width_ratios': [1, 3]})
    
    # Plot Sound Speed Profile
    if cProfile is None:
        ssp = [MunkSoundSpeedGrad(z, ssp_only=True) for z in sspDepths]
        fig.suptitle('Ocean Acoustic Ray Tracing with Munk Profile', fontsize=16)

    else:
        ssp = [SoundSpeedGrad(z, zProfile, cProfile,ssp_only=True) for z in sspDepths]
        fig.suptitle('Ocean Acoustic Ray Tracing with Given Profile', fontsize=16)

    ax1.plot(ssp, sspDepths)
    ax1.set_title("Sound Speed Profile")
    ax1.set_xlabel("Sound Speed (m/s)")
    ax1.set_ylabel("Depth (m)")
    ax1.invert_yaxis()
    ax1.grid(True)

    # Plot Ray Paths
    for angle, (path_data, travel_time) in ray_paths.items():
        ranges, depths = path_data
        label = f'{angle}°'
        if travel_time is not None:
             label += f' ({travel_time:.2f} s)'
        ax2.plot(ranges, depths, label=label)

    ax2.set_title("Acoustic Ray Paths")
    ax2.set_xlabel("Range (m)")
    ax2.set_ylabel("Depth (m)")
    ax2.invert_yaxis()
    ax2.legend(title="Initial Angle (Travel Time)")
    ax2.grid(True)
    #ax2.set_ylim(5000, -100) # Set depth limits for better visualization

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

def buildTTLookupTable(angleParams=[-90,-35,.5],depthParams=[0,6000,100],maxRange=50e3,zProfile=None,cProfile=None,dt=0.1, savetable=True):
    """
    Build a lookup table for travel times based on angle and depth.

    Args:
        angleParams (list): [min_angle, max_angle, angle_step]
        depthParams (list): [min_depth, max_depth, depth_step]
        maxRange (float): Maximum range for ray tracing (m).
        zProfile (array): Array of depths from the CTD profile (m).
        cProfile (array): Array of sound speeds corresponding to the depths (m/s).
        dt (float): Time step for the simulation (s).
        savetable (bool): Whether to save the lookup table to a file.

    Returns:
        dictionary: A dictionary mapping (range, initial_depth) to travel_time.
    """
    min_angle, max_angle, da = angleParams
    min_depth, max_depth, dz = depthParams
    if min_depth < dz:
        print(f"Adjusting min_depth from {min_depth} to {dz}")
        min_depth = dz
    # Define grid of angles and initial depths
    angles = np.arange(min_angle, max_angle + da, da)  # degrees
    source_depths = np.arange(min_depth, max_depth + dz, dz)
    # Define a function to run ray tracing
    def runRayTrace(args):
        source_depth, angle = args
        if zProfile is not None and cProfile is not None:
            ranges, depths, travel_time = RayTrace(angle, maxRange, source_depth, zProfile=zProfile, cProfile=cProfile, dt=dt)
        else:
            ranges, depths, travel_time = RayTrace(angle, maxRange, source_depth, dt=dt)
        if travel_time is not None:
            final_range = ranges[-1]
            return (final_range, source_depth, travel_time)
        else:
            return None
    # Lookup table: keys are (range, initial_depth), values are travel_time
    lookupxzt = {}
    print("Building lookup table for travel time by (range, initial_depth)... (parallelized)")
    # create a list of tasks to run over
    tasks = [(source_depth, angle) for source_depth in source_depths for angle in angles]
    try:
        from concurrent.futures import ThreadPoolExecutor, as_completed
        # Use ThreadPoolExecutor for numpy compatibility
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(runRayTrace, task) for task in tasks]
            for i, future in enumerate(as_completed(futures)):
                res = future.result()
                if res is not None:
                    final_range, source_depth, travel_time = res
                    lookupxzt[(final_range, source_depth)] = travel_time
                if (i+1) % 100 == 0 or (i+1) == len(tasks):
                    print(f"Progress: {i+1}/{len(tasks)}")
    except Exception as e:
        print(f"Parallel execution failed ({e}), falling back to serial.")
        for task in tasks:
            res = raytrace_task(task)
            if res is not None:
                final_range, source_depth, travel_time = res
                lookupxzt[(final_range, source_depth)] = travel_time
    if savetable:
        # Save the lookup table to a file
        with open("TTtable.pkl", "wb") as f:
            pickle.dump(lookupxzt, f)
    return lookupxzt

def interpXZT(query_range, query_depth, lookup_table, method='cubic'):
    """
    Interpolate to get travel time for given horizontal range and depth (both in meters).
    query_range, query_depth can be scalars or lists/arrays of same length.
    method: 'linear', 'nearest', or 'cubic' (see scipy.interpolate.griddata)
    Returns np.nan for points outside convex hull of the computed travel times.
    """
    lookup_points = np.array(list(lookup_table.keys()))  # shape (N, 2): columns are (range, depth)
    lookup_values = np.array(list(lookup_table.values()))  # shape (N,)
    # Convert inputs to arrays
    query_range = np.atleast_1d(query_range)
    query_depth = np.atleast_1d(query_depth)
    if query_range.shape != query_depth.shape:
        raise ValueError('query_range and query_depth must have the same shape')
    points = np.column_stack([query_range, query_depth])
    travel_time = griddata(lookup_points, lookup_values, points, method=method)
    # Print warning for any nan results
    for i, t in enumerate(travel_time):
        if np.isnan(t):
            print(f"Query point {points[i]} is outside the Travel Time lookup range")
    return travel_time if len(travel_time) > 1 else travel_time[0]