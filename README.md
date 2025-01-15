# A systematic approach for detecting abrupt shifts in ecological timeseries

This repository provides R scripts and data for the article: Pélissié M., Devictor V. & Dakos V. (2024) **A systematic approach for detecting abrupt shifts in ecological timeseries** published in *Biological Conservation* as part of the Special Issue *Non-equilibrium perspectives in biological conservation*. The article is available here: https://doi.org/10.1016/j.biocon.2023.110429.

## R

The R scripts have been implemented on R version 4.2.1.  
The 'R' directory contains functions that are used for the simulations and analyses.


## Data

The 'data' directory contains minimal data to generate the simulated timeseries and reproduce the analyses on ecological examples.


## Analyses

The R scripts are named according the order they should be sourced to reproduce the analyses and figures of the manuscript.  

`00_packages.R` lists and loads the packages used, if some that are not installed on your machine are not available on the CRAN the command to install them is indicated.  

`01_generate_timeseries.R` simulates and saves the timeseries to be used in analyses and figures.  

`02_run_classification.R` performs the trajectory classifications. **This step is long!** Optimally, it should be run on a multicore machine to parallel the computations.  

`03_figures.R` generates the material to make the figures from the manuscript's main text and supplementary.  


# How to classify your own timeseries?

A **graphical user interface** (GUI) is available [online](https://matpelissie.shinyapps.io/trajshift_app/) to classify single timeseries by providing a `csv` file.  
The `csv` file must consist in two columns, the first column indicates the time steps and the second column the variable of interest.  

The classification outputs highlight the best trajectory with the relevant metrics associated and can be compared with alternative trajectories. The quality of fit indices are indicated and compared to the ones from simulated timeseries (correctly and incorrectly classified).  

This app also allows to adjust the parameters used in the detection of breakpoints.
