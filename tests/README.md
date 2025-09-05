# Test Notebooks for OceanRays

The notebooks within show examples of running various codes within OceanRays

## Notebooks
1. `BuildSoundSpeed` shows how to:
    * Use existing CTD data in .cnv format to visualize variables
    * Build sound speed profiles from data
    * write out results to local file (and read them in again)
    * Compare results between different models and pre-computed values
 
1. `RayTracing` shows how to:
    * Run `RayTrace` to calculate the raypaths for acoustic waves through the above models from a starting depth and angle to either a maximum distance or a surface crossing.  The surface crossings are used later for travel-time calculations.  Examples are shown both for real data and a theoretical model.

1. `TTimes` show how to:
    * Use `buildTTLookupTable` to build a travel-time lookup-table for a specific model or CTD-derived velocity profile.
    * Using `interpXZT` to interrogate the built lookup-table and interpolate 1-way travel-times for a source/transponder depth and horizontal range.
    * Loading and using saved lookup-tables.
    * Contour plotting of the lookup-table data.