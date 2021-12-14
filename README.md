# expoDecayMod - R Package Readme

## Description
Package of implementation for exponential decay model in survival analysis

## Details
The Exponential Decay model allows the influence of an event to exponentially decay over time which is different from the original Cox model. In survival analysis, it can be used to investigate the trend of a treatment’s hazard that increases right after the treatment and then decreases until reaching the baseline hazard, such as organ transplant surgery.

This package contains a function for estimating the parameters in the Exponential Decay model.

## Usage
Install from GitHub, e.g., devtools::install_github("yiwen99/expoDecayMod")

Load the matrix package from the GUI or with the command library.

Type library(help="expoDecayMod") in R console to list the functions and data sets for the package. For more information on any of these, use ?(function_name). Use the example(function_name) for illustrations of usage.

## Reference
Keown-Stoneman, C., Horrocks, J., & Darlington, G. (2018). Exponential decay for binary time-varying covariates in Cox models. Statistics in medicine, 37(5), 776–788. 
