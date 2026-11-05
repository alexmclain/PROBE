# ------------------------------------------------------------------
# cnf_helpers.R
#
# Helper functions for the CNF (conditional normalizing flow) method
# used by spatial_lfdr_cnf() and e_step_func_sp() with method = "cnf".
#
# Functions:
#   lfdr_threshold()  — standard lfdr thresholding procedure
#   .cnf_build()      — internal: build MoL conditioner network
#   .mol_density()    — internal: evaluate mixture-of-logistics density
#   train_cnf()       — train a conditional normalizing flow (requires torch)
#   estimate_f_cnf()  — evaluate a fitted CNF at new (Z, coords) values
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# lfdr_threshold
#
# Standard local FDR thresholding: find the largest rejection set
# whose running average lfdr does not exceed alpha (Efron's procedure).
#
# Arguments:
#   lfdr_mat  : matrix or vector of estimated local FDR values
#   alpha     : target FDR level (default 0.1)
#   true_H    : optional true signal indicator (same shape as lfdr_mat)
#               for computing the empirical FDR
#   nrow_, ncol_ : dimensions for reshaping the rejection map; if NULL
#                  the rejection map is returned as a vector
#   verbose   : print number of discoveries (default TRUE)
#
# Returns a list:
#   alpha, k_star, R, V, empirical_FDR, estimated_FDR,
#   selected_indices, reject_mat, running_mean
# ------------------------------------------------------------------
lfdr_threshold <- function(lfdr_mat,
                           alpha    = 0.1,
                           true_H   = NULL,
                           nrow_    = NULL,
                           ncol_    = NULL,
                           verbose  = TRUE) {
  lfdr_vec <- as.numeric(lfdr_mat)
  n        <- length(lfdr_vec)
  
  ord          <- order(lfdr_vec)
  lfdr_sorted  <- lfdr_vec[ord]
  running_mean <- cumsum(lfdr_sorted) / seq_len(n)
  
  valid_k <- which(running_mean <= alpha)
  
  if (length(valid_k) == 0) {
    k_star   <- 0
    selected <- integer(0)
    if (verbose) cat("No discoveries at level alpha =", alpha, "\n")
  } else {
    k_star   <- max(valid_k)
    selected <- ord[seq_len(k_star)]
    if (verbose) cat("Number of discoveries:", k_star, "\n")
  }
  
  empirical_FDR <- NA
  V <- NA
  R <- k_star
  
  if (!is.null(true_H) && k_star > 0) {
    true_vec      <- as.numeric(true_H)
    V             <- sum(true_vec[selected] == 0)
    empirical_FDR <- V / R
    if (verbose) {
      cat("V (false discoveries):", V, "\n")
      cat("Empirical FDR:", empirical_FDR, "\n")
    }
  }
  
  reject_vec          <- rep(0L, n)
  reject_vec[selected] <- 1L
  
  if (!is.null(nrow_) && !is.null(ncol_)) {
    reject_mat <- matrix(reject_vec, nrow = nrow_, ncol = ncol_, byrow = FALSE)
  } else {
    reject_mat <- reject_vec
  }
  
  list(alpha          = alpha,
       k_star         = k_star,
       R              = R,
       V              = V,
       empirical_FDR  = empirical_FDR,
       estimated_FDR  = if (k_star > 0) running_mean[k_star] else 0,
       selected_indices = selected,
       reject_mat     = reject_mat,
       running_mean   = running_mean)
}


# ------------------------------------------------------------------
# .cnf_build  (internal)
#
# Build a fully-connected conditioner network:
#   input (2 spatial coords) -> hidden layers -> output (3K params)
# where K is the number of mixture-of-logistics components.
# Requires torch.
# ------------------------------------------------------------------
.cnf_build <- function(K, hidden_dim, n_hidden) {
  torch <- requireNamespace("torch", quietly = TRUE)
  if (!torch) stop("Package 'torch' is required. Install with install.packages('torch').")
  nn_linear   <- torch::nn_linear
  nn_tanh     <- torch::nn_tanh
  nn_sequential <- torch::nn_sequential
  
  layers <- list(nn_linear(2L, hidden_dim), nn_tanh())
  for (l in seq_len(n_hidden))
    layers <- c(layers, list(nn_linear(hidden_dim, hidden_dim), nn_tanh()))
  layers <- c(layers, list(nn_linear(hidden_dim, 3L * K)))
  do.call(nn_sequential, layers)
}


# ------------------------------------------------------------------
# .mol_density  (internal)
#
# Evaluate a K-component mixture-of-logistics density.
#   z      : (n, 1) torch tensor
#   params : (n, 3K) torch tensor — columns [log_pi | mu | log_s]
# Returns: (n,) torch tensor of density values (>= 0).
# ------------------------------------------------------------------
# .mol_density <- function(z, params, K) {
#   nnf_softmax  <- torch::nnf_softmax
#   nnf_softplus <- torch::nnf_softplus
#   
#   pi_k <- nnf_softmax(params[, 1:K],                        dim = 2L)
#   mu   <- params[, (K + 1):(2 * K)]
#   s_k  <- nnf_softplus(params[, (2 * K + 1):(3 * K)]) + 1e-3
#   
#   u   <- (z$expand(c(-1L, K)) - mu) / s_k
#   sig <- torch::torch_sigmoid(u)
#   f_k <- sig * (1 - sig) / s_k
#   
#   (pi_k * f_k)$sum(dim = 2L)
# }

.mol_density <- function(z, params, K) {     
  
  pi_k <- torch::nnf_softmax(params[, 1:K], dim = 2L)
  mu   <- params[, (K + 1):(2 * K)]
  s_k  <- torch::nnf_softplus(params[, (2 * K + 1):(3 * K)]) + 1e-3
  
  u   <- (z$expand(c(-1L, K)) - mu) / s_k   
  f_k <- torch::torch_exp(-0.5 * u^2) / (s_k * sqrt(2 * pi))                               
  
  (pi_k * f_k)$sum(dim = 2L)                    
}   


# ------------------------------------------------------------------
# train_cnf
#
# Train a conditional normalizing flow (mixture-of-logistics network)
# on all (Z_i, x_i, y_i) observations. The network maps spatial
# coordinates to the parameters of a K-component mixture of logistic
# densities: f(z | x,y) = sum_k pi_k(x,y) * Logistic(z; mu_k, s_k).
#
# Requires the 'torch' package.
#
# Arguments:
#   Z, coords   : nrow x ncol matrix and corresponding n x 2 coord matrix
#   K           : mixture components (default 10)
#   hidden_dim  : conditioner network width (default 64)
#   n_hidden    : conditioner hidden layers (default 2)
#   n_epochs    : Adam training epochs (default 1000)
#   lr          : Adam learning rate (default 1e-3)
#   batch_size  : mini-batch size (default 256)
#   verbose     : print progress every 500 epochs (default TRUE)
#
# Returns a list: net, K, cx_range, cy_range
# ------------------------------------------------------------------
train_cnf <- function(Z, coords,
                      K          = 10,
                      hidden_dim = 64,
                      n_hidden   = 2,
                      n_epochs   = 1000,
                      lr         = 1e-3,
                      batch_size = 256,
                      verbose    = TRUE) {
  if (!requireNamespace("torch", quietly = TRUE))
    stop("Package 'torch' is required. Install with install.packages('torch').")
  
  torch_tensor  <- torch::torch_tensor
  torch_float   <- torch::torch_float
  optim_adam    <- torch::optim_adam
  torch_log     <- torch::torch_log
  
  Zvec <- as.numeric(Z)
  n    <- length(Zvec)
  
  cx_range <- range(coords[, 1]); cy_range <- range(coords[, 2])
  cx_n <- 2 * (coords[, 1] - cx_range[1]) / diff(cx_range) - 1
  cy_n <- 2 * (coords[, 2] - cy_range[1]) / diff(cy_range) - 1
  
  z_ten    <- torch_tensor(Zvec,              dtype = torch_float())
  cond_ten <- torch_tensor(cbind(cx_n, cy_n), dtype = torch_float())
  
  net <- .cnf_build(K, hidden_dim, n_hidden)
  opt <- optim_adam(net$parameters, lr = lr)
  
  for (epoch in seq_len(n_epochs)) {
    idx     <- sample(n, min(batch_size, n))
    z_b     <- z_ten[idx]$unsqueeze(2L)
    c_b     <- cond_ten[idx, , drop = FALSE]
    params  <- net(c_b)
    density <- .mol_density(z_b, params, K)
    loss    <- -torch_log(density + 1e-12)$mean()
    
    opt$zero_grad()
    loss$backward()
    opt$step()
    
    if (verbose && epoch %% 500 == 0)
      cat(sprintf("CNF epoch %4d/%d  NLL = %.4f\n", epoch, n_epochs, loss$item()))
  }
  
  list(net = net, K = K, cx_range = cx_range, cy_range = cy_range)
}


# ------------------------------------------------------------------
# estimate_f_cnf
#
# Evaluate a trained CNF model at (Z, coords) and return a numeric
# vector of density values f(Z_i | x_i, y_i) >= 1e-12.
#
# Arguments:
#   Z, coords  : same or new data (coords must be within training range)
#   model      : object returned by train_cnf()
#   batch_size : evaluation batch size (default 512)
# ------------------------------------------------------------------
estimate_f_cnf <- function(Z, coords, model, batch_size = 512) {
  if (!requireNamespace("torch", quietly = TRUE))
    stop("Package 'torch' is required. Install with install.packages('torch').")
  
  torch_tensor     <- torch::torch_tensor
  torch_float      <- torch::torch_float
  with_no_grad     <- torch::with_no_grad
  
  Zvec <- as.numeric(Z)
  n    <- length(Zvec)
  K    <- model$K
  
  cx_n <- 2 * (coords[, 1] - model$cx_range[1]) / diff(model$cx_range) - 1
  cy_n <- 2 * (coords[, 2] - model$cy_range[1]) / diff(model$cy_range) - 1
  
  z_ten    <- torch_tensor(Zvec,              dtype = torch_float())
  cond_ten <- torch_tensor(cbind(cx_n, cy_n), dtype = torch_float())
  
  model$net$eval()
  f_hat_vec <- numeric(n)
  
  with_no_grad({
    for (start in seq(1, n, by = batch_size)) {
      idx              <- start:min(start + batch_size - 1L, n)
      params           <- model$net(cond_ten[idx, , drop = FALSE])
      density          <- .mol_density(z_ten[idx]$unsqueeze(2L), params, K)
      f_hat_vec[idx]   <- as.numeric(density)
    }
  })
  
  pmax(f_hat_vec, 1e-12)
}
