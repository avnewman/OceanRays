import numpy as np
from scipy.interpolate import griddata
import pickle

# local
from OceanRays.RayPaths import RayTrace

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
            res = runRayTrace(task)
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


## TODO
def readTTtablepkl():
    return None