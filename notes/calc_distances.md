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

1. Each __LC phase__  in each band have to be shifted to have its __zero-point__ 
   corresponding to the __phase at maximum flux in r-band__.

   This can be done before distance calculation. Call this **shifted-phase**.

2. When calculating the distance, **we have to know if the 2 LCs can be 
   compared**. This means that their relative shifted-phases have to overlap
   for some fraction.

3. This determination have to be done for each LC with respect to all the others.
   - The number of LC to compare to diminishes the deeper I go into the list of
     LC: 
     1. The first have to be compared to all the other
     2. The second to all but the first one 
     3. The third to all but the first and the second and so on
     4. **Can this be a recursive function?**

   - We'll end up with a NxN triangular matrix (upper or lower does not matter).

     Each cell will have a pair-wise distance and can be feeded to R's
     DiffusionMap package.


  