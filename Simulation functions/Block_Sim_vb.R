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

setwd("C:/Users/mclaina/OneDrive - University of South Carolina/Research/Imaging/PROBE/Programs/PROBE/Simulation functions")

remove(list=ls())
source("Sim_funcs.R")

# library(probe)
library(RcppArmadillo)
require(Rcpp)
require(ggplot2)
require(Matrix)
require(MASS)
require(glmnet)

library(dplyr)
library(doParallel)  

no_cores <- 28 # detectCores()
registerDoParallel(cores = no_cores)  
cl <- makeCluster(no_cores) 

parlist <- data.frame(read.csv("arg_list_block.csv",header = TRUE, fileEncoding="UTF-8-BOM"))
parlist <- parlist[order(-parlist$M),]
rownames(parlist) <- NULL

settings <- data.frame(seed = rep(seq(1000,10500,500),each=nrow(parlist)),parlist)
settings_to_try <- 1:nrow(settings)

foreach(set_i = settings_to_try)  %dopar% {
  
  library(glmnet)
  library(sparsevb)
  print(paste("Set =",set_i))
  args_list <- settings[set_i,-1]
  seed <- settings[set_i,1]
  output <- "C:/Users/mclaina/OneDrive - University of South Carolina/Research/Imaging/PROBE/Programs/PROBE_non_git/Sim Results/Results/Block/"
  sim1 <- block_simulation_VB(args_list, B = 20, ebreg_I = FALSE, verbose = TRUE, seed = seed, output = output)
  
  }
gc()
print("Finished.")

