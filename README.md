[![DOI](https://zenodo.org/badge/14986838.svg)](https://zenodo.org/badge/latestdoi/14986838)

miniature-adventure
===================

A photometric supernova classifier using a data driven approach.

Trained and tested on simulated dataset <http://sdssdp62.fnal.gov/sdsssn/SIMGEN_PUBLIC/SIMGEN_PUBLIC_DES.tar.gz>

Developed in Python and R languages as PhD project, and inspired to paper from [Richards, J.W. *et. al* (2012)](http://mnras.oxfordjournals.org/cgi/doi/10.1111/j.1365-2966.2011.19768.x). 

Thesis will be available soon after defense in fall 2015.


#### Processing steps:
1. Correction for astrophysical effects
2. Interpolation using Gaussian processes
3. Parameter extraction performed by diffusion maps
4. Classification model built with random forest algorithm

#### Important notes
- Interpolation is performed using Python package GPy, available at <https://github.com/SheffieldML/GPy>
- Parameter extraction is performed by R package diffusionMap
- The classification model is built using R package randomForests
