# OceanRays

## About

This github provides `OceanRays`; a set of python codes for the calculation of raypaths and travel times of acoustic waves in the open ocean.

The code will either use a simple theoretical ocean sound speed model, or can take in an external model, or create one from a CTD profile (using the .cnv format)

## Install 
### Create a usable Conda environment for working with the package
```
# create with a bunch of packages
conda create --name oceanrays python ipython ipykernel numpy scipy matplotlib
```
### Once you're in the proper environment, for me it is done by:
``` 
conda activate oceanrays
```
### Now install seabird for reading CTD data and tqdm for monitoring progress
```
pip install seabird tqdm
```
### Install for general user:
```
# from this directory
pip install .
# or from another directory
pip install '/path/to/this/directory/'
```
### Updating
```
pip install update . # from here
```
----
## Run
Look at the `.ipynb` notebooks within the `tests` directory for how to process. 