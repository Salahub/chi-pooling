# Controlling Centrality with the Chi-Squared Pooled p-Value

This repository contains all of the code required to run a recursive
binning algorithm to measure dependence. Two directories are contained:

## pooledCentR

The pooledCentR directory contains an R package which can be used to
estimate the centrality quotient, marginal rejection level at b, and
central rejection level for arbitrary p-value pooling functions. It
also contains an implementation of the chi-squared pooled p-value and
the uniformly most powerful (UMP) pooled p-value for a particular
restricted beta family. Supporting functions to generate examples
under different alternative distributions and visualize results are
provided alongside this core functionality. The package can be
downloaded via the command `install_github("Salahub/chi-pooling",
subdir="pooledCentR")` using `install_github` from the `devtools`
package in R.

## Scripts

This directory contains all extraneous data and scripts used to
undertake the explorations of centrality and the chi-squared pooled
p-value in my thesis. For tidiness, simulation outputs (which are
often called by other scripts) are stored in the results subdirectory.
Briefly:

- **poolFunctions.R**: essentially a prototype of pooledCentR, this
  contains similar functions to those found in the package along with
  many other helpful functions used throughout these scripts
- **pooledChi.R**: reads in many simulation results and generates
  plots to communicate these results, includes some rejection boundary
  and power curve investigations not included in my thesis
- **betaCentrality.R**: generates side-by-side plots of beta densities
  and the corresponding curves of the chi-square pooled p-value by
  kappa to demonstrate the relationship between the two (these plots
  are seen in Chapters 4 and 6) and simulate the null many times to
  generate reference quantiles
- **chiMapPlot.R**: processes power data from simulations using the
  chi-squared pooled p-value to create the data for functionality in
  pooledCentR
- **chiMap.R**: computes the power of the chi-square pooled p-value
  for a range of kappa values over a range of parameter settings (the
  simulation described in Section 4.7.2 of my thesis)
- **metaPooling**: investigates the use of the chi-squared pooling
  function to combine parameter estimates in meta-analysis through
  numerous simulation studies and a real data example
- **otherMethods**: implements several other pooled p-values along
  with functions that display their rejection regions in the case of
  two p-values being combined

The other script files (**simulateConstM_fullgrid.R**,
**simulateConstM.R**, **simulateVarM_fullgrid.R**, **simulateVarM.R**)
are jobs meant to be run on a machine without human intervention, and
so are generally not very easily read. They perform power simulations
from Sections 4.4 and 4.6.2 of my thesis which compare the
chi-squared pooled p-value to the uniformly most powerful pooled
p-value over a number of settings.
