import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import pickle

# pip install seabird (read CTD data from UNOLS)
from seabird.cnv import fCNV

# locals
from OceanRays.SoundSpeed import SoundSpeedGrad, MunkSoundSpeedGrad

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
# Plot sound speed profile and Ray Paths
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