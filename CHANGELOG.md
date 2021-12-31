## 31 December 2021 - v0.1.1
- Stopped using global namespace for several package variable. Defined
  Hickory environment internally.

## 26 December
- Added vignettes for locus- and population-specific effects and for
  dominant markers.

## 15 November
- Added outlier detection. NOTE: The approach used is different from
  the one described in Guo et al. It uses posterior comparison of
  individual effects to overall effects, taking advantage of the new
  parameterization.

## 8 November 2018
- Added estimates of locus- and population-specific effects on
  theta. NOTE: The parameterization is different from the one used in
  Guo et al. _Journal of the American Statistical Association_
  104:142-154; 2009.

## 18 October 2020

- Added estimate of within-population inbreeding coefficient by locus
  and population. NOTE: Each population locus combination is analyzed
  independently.
  
- Included summary of within-population inbreeding coefficient
  analysis to identify locus/population combinations in which the
  estimated credible interval (95% by default) for the
  within-population inbreeding coefficient does not include 0.
  
- Updated Roadmap. Defer investigation of dominant markers until after
  other tasks are complete.

## 11 October 2020

- Added vignette illustrating data format and analysis of data. 

- <img
   src="https://render.githubusercontent.com/render/math?math=\theta=0">
   is implemented.

## 3 October 2020

- <img src="https://render.githubusercontent.com/render/math?math=f=0">
   is implemented along with leave-one-out cross validation for model
   comparison.


   
