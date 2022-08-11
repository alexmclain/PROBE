### The following code can be use to reproduce the simulations given in the paper.
### The simulations were ran for all 54 combinations of parameters given in the 
### arg_list.csv file for binary and continuous Z. See the function "simulation_func" 
### in the "Sim_funcs.R" file for the code that runs the simulation.
### Note that the ebreg function takes markedly longer to run than the other 
### functions. To run the simulation without ebreg, use ebreg_I = FALSE when running
### simulation_func. In the paper, all methods were ran for 1000 iterations except
### ebreg, which was ran for 100 iterations. Also, ebreg was not ran for M=400 due 
### to the program frequently returning an error in those settings.
### simulation_func will automatically create a csv file with the results, but will 
### also return all of the results.
### The methods use are:
###   + UNHIDEM
###   + Penalization methods: LASSO, Adaptive LASSO, SCAD, MCP
###   + SSLASSO - Spike-and-Slab LASSO 
###   + varbvs -  Fully-factorized variational approximation for Bayesian variable 
###               selection
###   + sparsevb - Spike-and-Slab Variational Bayes
###   + sparsevb_c - Spike-and-Slab Variational Bayes with calibration (discussed in 
###                  paper)
###   + ebreg - Empirical Bayes sparse linear regression.


remove(list=ls())
source("Simulation functions/Sim_funcs.R")
library(unhidem)
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
colnames(parlist) <- c("sqrt M", "pi", "eta", "SNR", "Start")



### Run for continuous with M=400 for B=10 iterations without ebreg
par_num <- 5 
parlist[par_num,]
sim1 <- simulation_func(par_num, parlist, B = 10, bin = FALSE, ebreg_I = FALSE)

### Run for binary with M=2500 for B=10 iterations without ebreg
par_num <- 30 
parlist[par_num,]
sim2 <- simulation_func(par_num, parlist, B = 10, bin = TRUE, ebreg_I = FALSE)

### Run for continuous with M=10000 for B=10 iterations without ebreg.
### The difference between sparsevb and sparsevb_c is very evident here.
par_num <- 47 
parlist[par_num,]
sim3 <- simulation_func(par_num, parlist, B = 10, bin = FALSE, ebreg_I = FALSE)

### Run for binary with M=2500 for B=3 iterations with ebreg
par_num <- 28 
parlist[par_num,]
sim4 <- simulation_func(par_num, parlist, B = 3, bin = FALSE, ebreg_I = TRUE)
