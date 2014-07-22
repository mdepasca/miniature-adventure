miniature-adventure
===================

A photometric supernova classifier for SUDARE survey

- USES:
  - GPy devel branch (https://github.com/SheffieldML/GPy/tree/devel/GPy) 

- TODO:
  - Work on how to store on memory the fitting
  - KeyBoardInterrupt does not work
  - Investigate crash in parallel version of optimize_restart (worth it, parallel is 10 sec faster per 6 candidates)
	No more process deadlock in GPy/devel 12/07/2014 version.
  - How to use cleese for parallel computing.
