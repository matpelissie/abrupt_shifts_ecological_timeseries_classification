# A systematic approach for detecting abrupt shifts in ecological timeseries

This repository provides R scripts and data for the following manuscript: "A systematic approach for detecting abrupt shifts in ecological timeseries".

## R

The R scripts have been implemented on R version 4.2.1.  
The 'R' directory contains functions that are used for the simulations and analyses.


## Data

The 'data' directory contains minimal data to generate the simulated timeseries and reproduce the analyses on ecological examples.


## Analyses

The R scripts are named according the order they should be sourced to reproduce the analyses and figures of the manuscript.  
00_packages.R lists and loads the packages used, if some that are not installed on your machine are not available on the CRAN the command to install them is indicated.  
01_generate_timeseries.R simulates and saves the timeseries to be used in analyses and figures.  
02_run_classification.R performs the trajectory classifications. **This step is long!** Optimally, it should be run on a multicore machine to parallel the computations.  
03_figures.R generates the material to make the figures from the manuscript's main text and supplementary.  
