# Analysis of radiation damage to phospholipid monolayers

This repository provides the analysis method for a series of phospholipid monolayer (at the air-water interface) systems undergoing radiation damage.
The analysis performed herein has been developed by [Andrew R. McCluskey](mailto:andrew.mccluskey@diamond.ac.uk), with input from [Tom Arnold](mailto:tom.arnold@esss.se) and [Toby Robson](mailto:toby.robson@diamond.ac.uk).

This aim of this work is to investigate the mechanism of interaction between antioxidant additives and phospholipid monolayers.

## Model

The model for the phospholipid monolayer is as follows:

```
      air 
---------------
     tails
---------------
     heads
---------------
     water
```
Where there is interest in how the parameters change over time as they are exposed to the intense synchrotron radiation beam. 

## Bayesian Evidence

First, the Bayesian evidence for all permuations of the following parameters were found, tail thickness (`tt`), roughness (`rough`), tail solvation (`phit`) and head solvation (`phih`) was found.
These parameters were either constrainted to a single value as the radiation flux increases or allowed to vary. 
The following constraints were imposed on **all** analyses:
- Head thickness - `th = 10.` (this is constrained due to correlation between the head and tail thicknesses in XRR analysis)
- Head Volume - `mvh = 319.`
- Tail Volume - `mvt = 829.`
- Roughness was taken to be conformal. 
- Scale - `scale = 1.`
- Background - `bkg = dataset.y.min()` (minimum measured reflectometry value)

The Bayesian evidence was found as discussed in [Sivia and Skilling](https://global.oup.com/academic/product/data-analysis-9780198568322?cc=gb&lang=en&) for the monolayers 1, 4, 12, 22, and 33.
This is plotted below 

![]()

From this the ?? permutation was applied to all of the monolayers. 

## Time resolved analysis
