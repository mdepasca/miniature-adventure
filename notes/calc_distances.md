Calculating distances
======================

- To get the distance between 2 lightcurves is **easy**: R function `dist()`
  does excactly this

  From R help

  > Description:

  > This function computes and returns the distance matrix computed by
  > using the specified distance measure to compute the distances
  > between the rows of a data matrix.

Before calculating the distance
-------------------------------

- Each __LC phase__  in each band have to be shifted to have its __zero-point__ 
  corresponding to the __phase at maximum flux in r-band__
  