# CaveXC

A modular code for channel cross-section evolution in both free-surface and conduit-full (i.e. water filled caves) conditions.

This code provides a CrossSection class holding the geometry of the channel as a series of x, y points, with methods for calculation of area, wetted perimeter, position finding along the geometry, and evolution of the geometry. It also provides a method for finding a stage height given geometry and a discharge. A flow class is provided to calculate boundary shear stress along the cross-section based on the Wobus, Tucker, and Anderson (WTA) method from Wobus et al. [2006, 2008].

Two example scripts are included to simulate different cross-sections. The script vadose.py simulates surface bedrock channels or caves with open channel flow. This script requires a prescribed sediment diameter, slope, discharge, and initial geometry. The script paragenesis.py simulates paragenetic canyons, a unique type of cave passage where the channel grows upwards towards the water table in phreatic conditions when sediment armors the floor and walls of a passage from dissolution/mechanical erosion. As such, there is an additional module for calculating sediment transport dynamics, sediment.py.

The code paragenesis.py includes command line options given by "python paragenesis.py option=value". Physical parameters for the model are:

  * Q - discharge [m<sup>3</sup>/s],
  * Qs - sediment supply [kg/s],
  * rh - roughness length for the law of the wall (z<sub>0</sub>) [m],
  * D_s - sediment size [m],
  * n - exponent in the erosion law kτ<sub>b</sub><sup>n</sup>,
  * rho_s - sediment density [kg/m<sup>3</sup>].

Additionally, there are several option to control the model dynamics:
  * rh_type - type of roughness used in the model,
  * sc - whether suspended sediment is included in the model (0 - no suspended sediment/1 - suspended sediment),
  * t - number of time steps to run.

The option rh_type has several values:
  * 1 - use bedrock roughness as specified by physical parameter rh,
  * 2 - use sediment roughness calculated from D_s (6.8D_s/30),
  * 3 - composite roughness weighting bedrock roughness and sediment roughness by their perimeter portion lengths,
  * 4 - per-point roughness where z<sub>0</sub> varies during calculation of z<sub>0</sub>.

Further, options controlling the outputs of the script are:
  * sFile - file to save equilibrium cross-section data (eq. width, area),
  * suppress_paramsave - do not save equilibrium cross-section data (value=1),
  * save_tsdata - save per time step shear stress, cross-sectional area and width (value=1),
  * suppress_print - do not print per time step information to the command line.

Another script, runfromcsv.py, runs paragenesis.py given a CSV of options. Example CSVs are provided in examples/CooperCovington2020/. These CSVs are the options used for simulations in Cooper and Covington, in press, "Modeling cave cross-section including sediment transport and paragenesis", Earth Surface Processes and Landforms. Options for this script are:
  * sFile - the equivalent of sFile for paragenesis.py,
  * trf - CSV of options,
  * sedpart - the equivalent of sc for paragenesis.py.



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3401390.svg)](https://doi.org/10.5281/zenodo.3401390)
