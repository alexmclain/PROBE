### The following code can be use to reproduce the simulations given in the paper.
### The simulations were ran for all 54 combinations of parameters given in the 
### arg_list.csv file for binary and continuous X. See the function "simulation_func" 
### in the "Sim_funcs.R" file for the code that runs the simulation.
### Note that the ebreg function takes markedly longer to run than the other 
### functions. To run the simulation without ebreg, use ebreg_I = FALSE when running
### simulation_func. In the paper, all methods were ran for 1000 iterations except
### ebreg, which was ran for 100 iterations. Also, ebreg was not ran for M=400 due 
### to the program frequently returning an error in those settings.
### simulation_func will automatically create a csv file with the full results which
### will also be return. Additionally, the function will summarized results by
### iteration (if verbose=TRUE, the default) and for all iterations.
### The methods use are:
###   + PROBE (all-at-once and one-at-a-time)
###   + Penalization methods: 
###         -LASSO, 
###         -Adaptive LASSO, 
###         -SCAD, 
###         -MCP,
###   + SSLASSO - Spike-and-Slab LASSO 
###   + varbvs -  Fully-factorized variational approximation for Bayesian variable 
###               selection
###   + sparsevb - Spike-and-Slab Variational Bayes
###   + ebreg - Empirical Bayes sparse linear regression.
###
### Note: the following uses the geoR package to generate the data, while the
### simulations in the paper used the the RandomFields package. RandomFields
### is no longer available on CRAN and has become difficult to load. geoR has
### the same capabilities and is run with the same parameters, however, it
### it much much slower. For M=2500 generation of the data takes ~3 min. For
### M=10000 generation of the data takes multiple hours (with RandomFields
### it took minutes).

remove(list=ls())
source("Simulation functions/Sim_funcs.R")
library(probe)
library(glmnet)
library(sparsevb)
library(varbvs)
library(EMVS) 
### EMVS was removed from CRAN on Dec-2022 
### available at https://ftp.eenet.ee/pub/cran/web/packages/EMVS/index.html 
library(SSLASSO)
library(ebreg)
library(R.utils)
library(ncvreg)
library(dplyr)

parlist <- data.frame(read.csv("Simulation functions/arg_list.csv",header = TRUE, fileEncoding="UTF-8-BOM"))


### Run for binary predictors with M=400 for B=10 iterations without ebreg
bin = TRUE
par_num <- 13
args_list <- parlist[parlist$bin == bin,][par_num,]
sim1 <- simulation_func(args_list, B = 10, ebreg_I = FALSE, verbose = TRUE)

### Run for continuous predictors with M=2500 for B=5 iterations without ebreg
bin = FALSE
par_num <- 30 
parlist[parlist$bin == bin,][par_num,]
args_list <- parlist[parlist$bin == bin,][par_num,]
sim1 <- simulation_func(args_list, B = 5, ebreg_I = FALSE)

### Run for continuous predictors with M=2500 for B=3 iterations with ebreg
bin = TRUE
par_num <- 24 
parlist[parlist$bin == bin,][par_num,]
args_list <- parlist[parlist$bin == bin,][par_num,]
sim1 <- simulation_func(args_list, B = 3, ebreg_I = TRUE, verbose = TRUE)

### Run for binary predictors with M=10000 for B=5 iterations without ebreg.
bin = TRUE
par_num <- 48 
parlist[parlist$bin == bin,][par_num,]
args_list <- parlist[parlist$bin == bin,][par_num,]
sim1 <- simulation_func(args_list, B = 5, ebreg_I = FALSE, verbose = TRUE)
