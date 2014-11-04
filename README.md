miniature-adventure
===================

A photometric supernova classifier

- Simulated dataset SIMGEN_PUBLIC_DES.tar.gz at http://sdssdp62.fnal.gov/sdsssn/SIMGEN_PUBLIC/
  have to be unpacked in directory *train_data*

- USES:
  - GPy devel branch (https://github.com/SheffieldML/GPy/tree/devel/GPy) 

- TODO:
  - Work on how to store on memory the fitting: PANDAS? _its data structures are
    similar to R ones_
  - **`KeyBoardInterrupt` does not work**
  - Investigate crash in parallel version of optimize_restart 
    (worth it, parallel is 10 sec faster per 6 candidates)

    > No more process deadlock in GPy/devel 12/07/2014 version.
  - How to use cleese for parallel computing.
