# OceanRays

## About

This github provides `OceanRays`; a set of python codes for the calculation of raypaths and travel times of acoustic waves in the open ocean.

The code will either use a simple theoretical ocean sound speed model, or can take in an external model, or create one from a CTD profile (using the .cnv format)

## Install 
### Create a usable Conda environment for working with the package
```
# create with a bunch of packages
conda create --name oceanrays python ipython ipykernel numpy scipy pickle
```
### Once you're in the proper environment, for me it is done by:
``` 
conda activate oceanrays
```
### Now install seabird for reading CTD data
```
pip install seabird
```

### Install the `OceanRays` package by going into the directory with `setup.py` and using pip:

#### Install for general user:
```
pip install .
# alternatively, install from another directory, but point to that directory
pip install '/path/to/setup/.py/file'
```
To update your install, after downloading (e.g. through git pull)
```
# remember to be in the `setup.py` directory
pip install update .
```
#### Install for developer:
Alternatively, if you are planning on modifying the code for your own development, it would likely be better to run the below, to allow updates to code to be automatically loaded

Run in 'Edit' Mode
```
# Allowing edits to programs to automatically become loaded
pip install -e .
```
----
## Run
Look at the `.ipynb` notebooks for how to process. 