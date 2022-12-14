### The following code can be use to reproduce the simulations given in the paper.
### The simulations were ran for all 54 combinations of parameters given in the 
### arg_list.csv file for binary and continuous Z. See the function "simulation_func" 
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
###   + PROBE
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


remove(list=ls())
source("Simulation functions/Sim_funcs.R")
library(probe)
library(glmnet)
library(sparsevb)
library(varbvs)
library(EMVS)
library(SSLASSO)
library(ebreg)
library(R.utils)
library(ncvreg)
library(dplyr)

parlist <- as.matrix(data.frame(read.csv("Simulation functions/arg_list.csv",header = FALSE, fileEncoding="UTF-8-BOM")))
colnames(parlist) <- c("sqrt M", "pi", "eta", "SNR")

 

### Run for binary predictors with M=400 for B=10 iterations without ebreg
par_num <- 13
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- simulation_func(args_list, B = 10, bin = TRUE, ebreg_I = FALSE)

### Run for continuous predictors with M=2500 for B=5 iterations without ebreg
par_num <- 30 
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- simulation_func(args_list, B = 5, bin = FALSE, ebreg_I = FALSE)

### Run for continuous predictors with M=10000 for B=5 iterations without ebreg.
par_num <- 48 
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- simulation_func(args_list, B = 5, bin = FALSE, ebreg_I = FALSE)

### Run for binary predictors with M=2500 for B=3 iterations with ebreg
par_num <- 23 
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- simulation_func(args_list, B = 3, bin = TRUE, ebreg_I = TRUE)
