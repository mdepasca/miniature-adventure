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


Looking for overlapping regions
-------------------------------
- I keep light curves sampled on a 1-MJD grid; every sampling **is not** at a 
  a integer value of an MJD.

  When it will be time to compare fluxes, they will be probably associated to 
  different MJDs, but the difference will be smaller then 1 MJD and this 
  is most likely not to produce an significant error. **Anyway, this could be
  teste**.

- **The overlapping region**: if the maximums of the 2 LC are contained in the 
  respective phase limits. There will be 2 overlapping arrays, one per LC.
    - *Their number of elements*: will be the sum of 
      1. the number of elements before the zero-phase happening first
      2. the number of elements after the zero-phase happening last

      In each of the 2 the zero-phase point is included (it is the MJD 
      associated with the maximum in r-band).

1. get the elements of the first array that are in the second
2. get the elements of the second array that are in the first

   These two instruction are to build 2 masks. They are performed using
   `np.in1d()` function with `invert=True` option (so that to have a mask)
3. As per construction, the unmasked values will be matching MJDs.
4. The mask is applied to flux and flux errors arrays and the differences
   can be performed 
  