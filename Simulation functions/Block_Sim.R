### The following code can be use to reproduce the simulations given in the paper.
### The simulations were ran for all 30 combinations of parameters given in the 
### arg_list_block.csv file for binary and continuous X. See the function "block_simulation" 
### in the "Sim_funcs.R" file for the code that runs the simulation.
###
### simulation_func will automatically create a csv file with the full results.
### Additionally, the function will summarized results by iteration (if 
### verbose=TRUE, the default) and for all iterations.
###
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


remove(list=ls())
source("Sim_funcs.R")
library(probe)
library(glmnet)
library(sparsevb)
library(EMVS) 
### EMVS was removed from CRAN on Dec-2022 
### available at https://ftp.eenet.ee/pub/cran/web/packages/EMVS/index.html 
library(SSLASSO)
library(ebreg)
library(R.utils)
library(ncvreg)
library(dplyr)

## File path to output simulation results
output <- NULL

parlist <- data.frame(read.csv("arg_list_block.csv",header = TRUE, fileEncoding="UTF-8-BOM"))


par_num <- 1
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- block_simulation(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)


par_num <- 15 
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- block_simulation(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)


par_num <- 30 
parlist[par_num,]
args_list <- parlist[par_num,]
sim1 <- block_simulation(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)



### Robust Simulations


par_num <- 1
parlist[par_num,]
## Add degrees of freedom for the t-distribution
args_list <- c(parlist[par_num,],10) 
sim1 <- block_simulation_robust(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)


par_num <- 15 
parlist[par_num,]
args_list <- parlist[par_num,]
## Add degrees of freedom for the t-distribution
args_list <- c(parlist[par_num,],10) 
sim1 <- block_simulation_robust(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)


par_num <- 30 
parlist[par_num,]
args_list <- parlist[par_num,]
## Add degrees of freedom for the t-distribution
args_list <- c(parlist[par_num,],10) 
sim1 <- block_simulation_robust(args_list, B = 3, ebreg_I = FALSE, verbose = TRUE, output = output)



