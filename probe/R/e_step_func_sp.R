#' @title Spatial empirical Bayes E-step using spatially-weighted KDE and a local truncated normal null
#' @description Estimates posterior expectations of the \eqn{\gamma} variables for the spatial PROBE
#' algorithm. For each predictor location, test statistics are pooled from spatial neighbors using
#' Gaussian kernel weights. The marginal density is estimated via a weighted KDE, and the null density
#' is estimated via a locally-weighted truncated normal model with variance constrained to be at least one.
#' Local \eqn{\pi_0} estimates follow the weighted Storey method.
#' @usage e_step_func_sp(beta_t, beta_var, df, Dfull_grid, neighbor_list,
#' h_spatial = NULL, h_z = NULL, radius = NULL, null_frac = 0.70,
#' lambda = 0.5, min_neighbors = 10, lambda_tau = 0.1, H = NULL, coord_leng)
#' @param beta_t Numeric vector of posterior mean expectations (assuming \eqn{\gamma=1}).
#' @param beta_var Numeric vector of current posterior variances (assuming \eqn{\gamma=1}).
#' @param df Degrees of freedom for the t-distribution used to compute p-values.
#' @param Dfull_grid An \code{M x M} matrix of pairwise distances between all predictor locations
#' (e.g., from \code{fields::rdist(coords)}).
#' @param neighbor_list A precomputed list of length \code{M}, where element \code{i} contains the
#' integer indices of spatial neighbors for location \code{i}.
#' @param h_spatial Spatial Gaussian kernel bandwidth. If \code{NULL}, uses a heuristic based on
#' the median pairwise distance.
#' @param h_z KDE bandwidth in the test statistic space. If \code{NULL}, uses a weighted Silverman
#' rule of thumb computed within each neighborhood.
#' @param radius Neighbor radius used to construct \code{neighbor_list} (informational; passed through
#' to diagnostic output).
#' @param null_frac Fraction of central (near-zero) test statistics used to define the local null
#' region (default = 0.70).
#' @param lambda Storey tuning parameter for local \eqn{\pi_0} estimation (default = 0.5).
#' @param min_neighbors Minimum number of neighbors required for local estimation. Locations with fewer
#' neighbors fall back to using all \code{M} predictors (default = 10).
#' @param lambda_tau Ridge-type penalty on the local null variance estimate to prevent instability
#' in sparse neighborhoods (default = 0.1).
#' @param H (optional) True signal matrix, used only for diagnostic image plots when \code{verbose = TRUE}
#' inside the inner function.
#' @param coord_leng Integer vector of coordinate dimension lengths (e.g., \code{c(nx, ny)} for a
#' 2D grid) used to reshape results into matrices for verbose diagnostic plots.
#' @return A list including
#'
#' \code{gamma} estimated posterior expectations of the \eqn{\gamma} variables,
#'
#' \code{lfdr} estimated local false discovery rates (\eqn{1 - \gamma}),
#'
#' \code{pi0} local estimates of the proportion of null hypotheses,
#'
#' \code{p_vals} two-sided p-values computed from the t-distribution,
#'
#' \code{T_vals} test statistics (\eqn{\hat{\beta} / \sqrt{\widehat{\mathrm{Var}}(\hat{\beta})}}),
#'
#' \code{beta_tilde} posterior mean expectations (equal to input \code{beta_t}),
#'
#' \code{beta_tilde_var} posterior variances (equal to input \code{beta_var}),
#'
#' \code{adj_warning} integer flag (1 if the marginal density appears too narrow relative to the null density, 0 otherwise).
#' @references 
#' \itemize{  \item McLain, AC, A Zgodic, H Bondell (2025). Sparse high-dimensional linear regression with a partitioned empirical Bayes ECM algorithm. \emph{Computational Statistics and Data Analysis} 207, 108146.
#' \item Zgodic, A., Bai, R., Zhang, J., Wang, Y., Rorden, C., & McLain, AC. (2026). Quantifying predictive uncertainty of aphasia severity in stroke patients with sparse heteroscedastic Bayesian high-dimensional regression. \emph{Computational Statistics} 41:7, 1 -- 23.
#' \item Storey, J. D., Taylor, J. E., and Siegmund, D. (2004), “Strong control, conservative point estimation and simultaneous conservative consistency of false discovery rates: A unified approach,” \emph{J. R. Stat. Soc. Ser. B. Stat. Methodol.}, 66, 187–205.
#' }
#' @examples
#' #not run
#' #mod <- e_step_func(beta_t, beta_var, df, adj = 5, lambda = 0.1, monotone = TRUE)

e_step_func_sp <- function(beta_t, beta_var, df,
                           Dfull_grid,
                           neighbor_list,
                           h_spatial = NULL,    # spatial kernel bandwidth (if NULL use heuristic)
                           h_z = NULL,          # KDE bandwidth in z (if NULL use Silverman-ish)
                           radius = NULL,       # neighbor radius in grid units (if NULL use 6)
                           null_frac = 0.70,    # fraction of central Z to treat as local null
                           lambda = 0.5,        # Storey tuning
                           min_neighbors = 10,  # fallback if too few neighbors
                           lambda_tau = 0.1,    # Penalty for estimate of sigma
                           H = NULL,
                           coord_leng
                           ){
  
  ### Get and clean up T statistics
  T_vals <- beta_t/sqrt(beta_var)
  T_vals[beta_var<=0] <- 0
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0   
  p_vals <- pt(abs(T_vals),df=df,lower.tail = FALSE)*2
  
  ### Estimate the lfdr
  EB_res <- spatial_nonparam_kde_truncnorm(T_vals,
                                           p_vals,
                                           Dfull_grid,
                                           neighbor_list, # list of neighbors
                                           h_spatial, 
                                           h_z,          # KDE bandwidth in z (if NULL use Silverman-ish)
                                           radius,       # neighbor radius in grid units (if NULL use 6)
                                           null_frac,    # fraction of central Z to treat as local null 
                                           lambda,        # Storey tuning
                                           min_neighbors,  # fallback if too few neighbors
                                           lambda_tau,     # Penalty for estimate of sigma
                                           df,
                                           H,            # True value of the signals
                                           eps = 1e-12,
                                           coord_leng,
                                           verbose = FALSE)
  gamma <- 1-EB_res$lfdr
  adj_warning <- EB_res$adj_warning

  
  ret <- list(beta_tilde = beta_t, beta_tilde_var = beta_var, gamma = gamma, 
              lfdr = 1 - gamma, pi0 = EB_res$pi0, p_vals = p_vals, T_vals = T_vals,
              adj_warning = adj_warning)
  return(ret)
}





spatial_nonparam_kde_truncnorm <- function(T_vals,
                                           p_vals,
                                           Dfull_grid,
                                           neighbor_list, # list of neighbors
                                           h_spatial = 0.1, 
                                           h_z = NULL,          # KDE bandwidth in z (if NULL use Silverman-ish)
                                           radius = NULL,       # neighbor radius in grid units (if NULL use 6)
                                           null_frac = 0.70,    # fraction of central Z to treat as local null 
                                           lambda = 0.5,        # Storey tuning
                                           min_neighbors = 10,  # fallback if too few neighbors
                                           lambda_tau = 0.1,     # Penalty for estimate of sigma
                                           df = 40,
                                           H = NULL,            # True value of the signals
                                           eps = 1e-12,
                                           coord_leng,
                                           verbose = TRUE) {
  
  n <- length(T_vals)
  
  # precompute p-values and global central threshold 
  z_thres <- as.numeric(stats::quantile(abs(T_vals), probs = null_frac, na.rm = TRUE))
  z_abs_q <- abs(T_vals) <= z_thres
  
  # precompute a global marginal KDE for fallback
  temp_h_z <- h_z
  if (is.null(h_z)) {
    # Silverman-ish bandwidth on T_vals (local)
    temp_h_z <- weighted_h(T_vals) 
  }
  global_kde <- stats::density(T_vals, bw = temp_h_z)
  global_fun <- approxfun(global_kde$x, global_kde$y, rule = 2)
  
  # allocate result vectors
  f_hat_vec <- numeric(n)
  f0_hat_vec <- numeric(n)
  pi0_vec <- numeric(n)
  lfdr_vec <- numeric(n)
  h_z_vec <- numeric(n)
  sigma_vec <- numeric(n)
  
  global_sigma <- wt_truncnorm_sigma_ge1(T_vals[z_abs_q], rep(1,sum(z_abs_q)), t = z_thres)
  if (verbose) {
    message("spatial_nonparam_kde: n = ", n, ", radius = ", radius,
            ", h_spatial = ", round(h_spatial,3), ", global null sigma = ", round(global_sigma, 3))
  }
  # loop locations (can be parallelized)
  for (i in seq_len(n)) {
    
    neigh <- neighbor_list[[i]]
    if (length(neigh) < min_neighbors) {
      # expand to full grid if too few neighbors
      neigh <- seq_len(n)
    }
    # spatial weights: Gaussian on distance^2
    d2 <- Dfull_grid[neigh,i]^2
    w_spatial <- exp(- d2 / (2 * h_spatial^2))
    if (sum(w_spatial) <= 0) w_spatial <- rep(1, length(w_spatial))
    # ---- marginal f: weighted KDE on observed Z in neighborhood ----
    samples_obs <- T_vals[neigh]
    
    temp_h_z <- h_z
    if (is.null(h_z)) {
      # Silverman-ish bandwidth on T_vals (local)
      temp_h_z <- weighted_h(samples_obs, w_spatial) 
    }
    f_hz <- temp_h_z
    f_hat_i <- tryCatch(weighted_kde_eval(samples_obs, w_spatial, T_vals[i], temp_h_z),
                        error = function(e) global_fun(T_vals[i]))
    # fallback
    if (!is.finite(f_hat_i) || f_hat_i <= 0) f_hat_i <- max(global_fun(T_vals[i]), eps)
    f_hat_vec[i] <- f_hat_i
    # ---- null f0 estimation ----
    # central mask among neigh
    central_mask <- z_abs_q[neigh]
    if (sum(central_mask) < min_neighbors) {
      # if too few central points, fall back to global central threshold across all neighbors, else to global KDE
      central_mask <- rep(TRUE, length(neigh))
    }
    samples_null <- T_vals[neigh][central_mask]
    weights_null <- w_spatial[central_mask]
    
    # Estimate sigma
    est_sigma <- wt_truncnorm_sigma_ge1(samples_null, weights_null, t = z_thres, lambda_tau = lambda_tau)
    
    if (length(samples_null) < 3) {
      # fallback: use global estimate of sigma
      f0_i <- dt(T_vals[i]/global_sigma, df = df)
      sigma_vec[i] <- global_sigma
    } else {
      f0_i <- dt(T_vals[i]/est_sigma, df = df)
      sigma_vec[i] <- est_sigma
    }
    
    if (!is.finite(f0_i) || f0_i <= 0) f0_i <- eps
    f0_hat_vec[i] <- f0_i
    # ---- local pi0 via weighted Storey ----
    p_neigh <- p_vals[neigh]
    w_p <- w_spatial
    numer <- sum(w_p * (p_neigh > lambda))
    denom <- (1 - lambda) * sum(w_p)
    pi0_i <- min(1, numer / max(denom, eps))
    pi0_vec[i] <- pi0_i
    # ---- lfdr plug-in and clamp ----
    lfdr_i <- min(1, pi0_i * f0_i / f_hat_i)
    lfdr_vec[i] <- max(0, min(1, lfdr_i))
    h_z_vec[i] <- f_hz
  } # end loop
  
  if(length(coord_leng) == 2){
    # reshape to matrices
    if (verbose) {
      lfdr_mat <- matrix(lfdr_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      f_mat <- matrix(f_hat_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      f0_mat <- matrix(f0_hat_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      pi0_mat <- matrix(pi0_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      h_z_mat <- matrix(h_z_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      sigma_mat <- matrix(sigma_vec, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      T_vals_mat <- matrix(T_vals, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
      op <- par(no.readonly = TRUE)
      par(mfrow = c(3,3), mar = c(2,2,2,1))
      if (!is.null(H)){
        fields::image.plot(t(H)[, nrow(H):1], main = "True H (mask)")
      }else{
        hist(p_vals, main = "Observed p histogram", breaks = 40, col = gray(0.9))
      }
      fields::image.plot(t(T_vals_mat)[, nrow(T_vals_mat):1], main = "Observed T")
      fields::image.plot(t(lfdr_mat)[, nrow(lfdr_mat):1], main = "Estimated lfdr")
      fields::image.plot(t(f0_mat)[, nrow(lfdr_mat):1], main = "Estimated f0")
      fields::image.plot(t(f_mat)[, nrow(lfdr_mat):1], main = "Estimated f")
      fields::image.plot(t(pi0_mat)[, nrow(lfdr_mat):1], main = "Estimated pi0")
      fields::image.plot(t(h_z_mat)[, nrow(lfdr_mat):1], main = "Bandwidth")
      fields::image.plot(t(sigma_mat)[, nrow(lfdr_mat):1], main = "Estimated Null Sigma")
      plot(sort(array(lfdr_mat)), type = "l", main = "Sorted lfdr values")
      par(op)
    }
  }
  
  ## Bandwidth warning:
  max_d <- global_fun(0)
  max_n <- c(f0_hat_vec, dt(0, df = df))
  adj_warning = 0
  if(any(max_d > max_n*1.25)){
    adj_warning = 1
  }
  
  
  return(list(lfdr = lfdr_vec,
              f = f_hat_vec,
              f0 = f0_hat_vec,
              pi0 = pi0_vec,
              h_z = h_z_vec,
              sigma = sigma_vec,
              params = list(h_spatial = h_spatial, h_z = h_z, radius = radius, null_frac = null_frac),
              neighbor_list = neighbor_list, 
              adj_warning = adj_warning))
}




# This uses a parametric normal local null with mean fixed at zero and 
# variance constrained to be at least one, estimated locally via weighted 
# likelihood, with shrinkage to stabilize the fit.
wt_truncnorm_sigma_ge1 <- function(z0, w0, t = 2, lambda_tau = 0.1) {
  
  neglog <- function(theta) {
    mu <- 0; sigma <- theta[1]
    a <- (-t - mu) / sigma; b <- (t - mu) / sigma
    Znorm <- pnorm(b) - pnorm(a)
    if (Znorm <= 1e-12) return(1e20)
    u <- (z0 - mu) / sigma
    logpdf <- dnorm(u, log = TRUE) - log(sigma) - log(Znorm)
    nll <- -sum(w0 * logpdf) + 0.5 * lambda_tau * (sigma - 1)^2
    nll
  }
  opt <- optimize(neglog, c(1, 20))
  sigma_hat <- opt$minimum
  sigma_hat
}


# weighted KDE evaluator (Gaussian kernel). samples: numeric vector, weights: same length
weighted_kde_eval <- function(samples, weights, x0, h, eps = 1e-12) {
  # small defensive checks
  if (length(samples) == 0) return(eps)
  # avoid zero total weight
  wsum <- sum(weights)
  if (wsum <= 0) return(eps)
  u <- (x0 - samples) / h
  # Gaussian kernel density
  ker_vals <- dnorm(u)
  val <- sum(weights * ker_vals) / (h * wsum)
  # guard
  if (!is.finite(val) || val <= 0) return(eps)
  return(val)
}

# Calculates the bandwidth using the Silverman rule with weights.
# The weights are used to calculate a weighted sd and the effective
# sample size.
weighted_h <- function(samples_obs, w = NULL) {
  if(is.null(w)) w <- rep(1,length(samples_obs))
  # Silverman-ish bandwidth on T_vals (local)
  sd_z <- max(c(weighted_sd(samples_obs, w),1))
  n_eff <- sum(w)
  temp_h_z <- 1.06 * sd_z * n_eff^(-1/5)
  if (temp_h_z <= 0 || is.na(temp_h_z)) temp_h_z <- 0.3
  temp_h_z
}

# Calculates a weighted standard deviation.
weighted_sd <- function(x, w) {
  mu <- sum(w * x) / sum(w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}

