#' @title Spatial empirical Bayes E-step using spatially-weighted KDE or CNF and a local truncated normal null
#' @description Estimates posterior expectations of the \eqn{\gamma} variables for the spatial PROBE
#' algorithm. For each predictor location, test statistics are pooled from spatial neighbors using
#' Gaussian kernel weights. The marginal density is estimated via either a weighted KDE
#' (\code{method = "kde"}, default) or a conditional normalizing flow neural network
#' (\code{method = "cnf"}), and the null density is estimated via a locally-weighted truncated normal
#' model with variance constrained to be at least one. Local \eqn{\pi_0} estimates follow the
#' weighted Storey method.
#' @usage e_step_func_sp(beta_t, beta_var, df, Dfull_grid, neighbor_list,
#' h_spatial = NULL, h_z = NULL, radius = NULL, null_frac = 0.70,
#' lambda = 0.5, min_neighbors = 10, lambda_tau = 0.1, H = NULL, coord_leng,
#' method = c("kde", "cnf"), coords = NULL, estimate_sigma = TRUE, sigma_null = 1,
#' cnf_model = NULL, cnf_K = 10, cnf_hidden_dim = 64, cnf_n_hidden = 2,
#' cnf_n_epochs = 1000, cnf_lr = 1e-3, cnf_batch_size = 256, sp_weights = NULL,
#' monotone = FALSE, monotone_step = 0.1)
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
#' @param method Density estimation method for the marginal distribution. \code{"kde"} (default)
#' uses spatially-weighted kernel density estimation. \code{"cnf"} uses a conditional normalizing
#' flow (neural-network mixture-of-logistics) and requires the \code{torch} package along with
#' the helper functions \code{train_cnf}, \code{estimate_f_cnf}, and \code{lfdr_threshold}.
#' @param coords (required when \code{method = "cnf"}) An \code{M x d} matrix of spatial
#' coordinates for each predictor location.
#' @param estimate_sigma Logical; only used when \code{method = "cnf"}. If \code{TRUE} (default),
#' estimates the local null variance via a weighted truncated-normal MLE with shrinkage. If
#' \code{FALSE}, fixes null SD to \code{sigma_null}.
#' @param sigma_null Fixed null SD used when \code{method = "cnf"} and \code{estimate_sigma = FALSE}
#' (default = 1).
#' @param cnf_model A previously fitted CNF model (the \code{$cnf_model} element from a prior call
#' with \code{method = "cnf"}). When supplied, training is skipped. Ignored for \code{method = "kde"}.
#' @param cnf_K Number of mixture-of-logistics components for the CNF (default = 10).
#' @param cnf_hidden_dim Conditioner network width for the CNF (default = 64).
#' @param cnf_n_hidden Number of hidden layers in the CNF conditioner (default = 2).
#' @param cnf_n_epochs Adam training epochs for the CNF (default = 1000).
#' @param cnf_lr Adam learning rate for the CNF (default = 1e-3).
#' @param cnf_batch_size Mini-batch size for CNF training (default = 256).
#' @param sp_weights (optional) Precomputed spatial weights from \code{build_sp_weights()}; only used
#' when \code{method = "kde"}. \code{sp_probe} computes these once per fit. If \code{NULL}, they are
#' computed from \code{neighbor_list}, \code{Dfull_grid} and \code{h_spatial}.
#' @param monotone Logical; only used when \code{method = "kde"}. If \code{TRUE}, each location's
#' lfdr is made non-increasing in \eqn{|T|}, as in \code{e_step_func(monotone = TRUE)}: the lfdr at
#' \eqn{T_j} is replaced by the maximum of location \eqn{j}'s lfdr curve over
#' \eqn{|T_j| \le |t| \le t_{max,j}} on the same side of zero as \eqn{T_j}, where \eqn{t_{max,j}} is
#' the largest \eqn{|T|} among \eqn{j}'s neighbors (so the tails of the local density, beyond the
#' local data, have no effect). Since the curve differs by location, it is evaluated on a grid with
#' spacing \code{monotone_step} (default = \code{FALSE}). This is considerably slower.
#' @param monotone_step Grid spacing in \eqn{t} used when \code{monotone = TRUE} (default = 0.1).
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
#'
#' @importFrom stats approxfun dnorm optimize pnorm
#' @importFrom graphics par

e_step_func_sp <- function(beta_t, beta_var, df,
                           Dfull_grid,
                           neighbor_list,
                           h_spatial = NULL,       # spatial kernel bandwidth (if NULL use heuristic)
                           h_z = 1,               # KDE bandwidth in z multiplier (if 1 uses Silverman-ish)
                           radius = NULL,          # neighbor radius in grid units (if NULL use heuristic)
                           null_frac = 0.70,       # fraction of central Z to treat as local null
                           lambda = 0.5,           # Storey tuning
                           min_neighbors = 10,     # fallback if too few neighbors
                           lambda_tau = 0.1,       # Penalty for estimate of sigma
                           H = NULL,
                           coord_leng,
                           verbose = FALSE,
                           method = c("kde", "cnf"),
                           coords = NULL,          # required for method = "cnf"
                           estimate_sigma = TRUE,  # CNF: estimate local null sigma (TRUE) or fix to sigma_null
                           sigma_null = 1,         # CNF: fixed null SD when estimate_sigma = FALSE
                           cnf_model = NULL,       # CNF: pre-fitted model to skip retraining
                           cnf_K = 10,
                           cnf_hidden_dim = 64,
                           cnf_n_hidden = 2,
                           cnf_n_epochs = 1000,
                           cnf_lr = 1e-3,
                           cnf_batch_size = 256,
                           sp_weights = NULL,      # KDE: precomputed spatial weights
                           monotone = FALSE,       # KDE: make each location's lfdr non-increasing in |T|
                           monotone_step = 0.1     # KDE: grid spacing for monotone
                           ){
  method <- match.arg(method)

  ### Get and clean up T statistics
  T_vals <- beta_t/sqrt(beta_var)
  T_vals[beta_var<=0] <- 0
  T_vals[is.na(T_vals)] <- 0
  T_vals[is.nan(T_vals)] <- 0
  p_vals <- pt(abs(T_vals), df = df, lower.tail = FALSE) * 2

  # ---- CNF branch ----
  if (method == "cnf") {
    if (is.null(coords))
      stop("coords must be provided when method = 'cnf'.")
    if (!requireNamespace("torch", quietly = TRUE))
      stop("Package 'torch' is required for method = 'cnf'. Install with install.packages('torch').")
    # Reshape T_vals to matrix expected by spatial_lfdr_cnf
    if (length(coord_leng) == 2) {
      T_mat <- matrix(T_vals, nrow = coord_leng[1], ncol = coord_leng[2])
    } else {
      T_mat <- matrix(T_vals, ncol = 1)
    }
    cnf_res <- spatial_lfdr_cnf(
      Z              = T_mat,
      coords         = coords,
      h_spatial      = h_spatial,
      Dfull_grid     = Dfull_grid,
      radius         = radius,
      null_frac      = null_frac,
      lambda         = lambda,
      min_neighbors  = min_neighbors,
      estimate_sigma = estimate_sigma,
      sigma_null     = sigma_null,
      lambda_tau     = lambda_tau,
      H              = H,
      eps            = 1e-12,
      cnf_model      = cnf_model,
      cnf_K          = cnf_K,
      cnf_hidden_dim = cnf_hidden_dim,
      cnf_n_hidden   = cnf_n_hidden,
      cnf_n_epochs   = cnf_n_epochs,
      cnf_lr         = cnf_lr,
      cnf_batch_size = cnf_batch_size,
      verbose        = verbose
    )
    lfdr_vec <- as.numeric(cnf_res$lfdr)
    gamma    <- 1 - lfdr_vec

    ret <- list(
      gamma       = gamma,
      lfdr        = lfdr_vec,
      pi0         = as.numeric(cnf_res$pi0),
      p_vals      = p_vals,
      T_vals      = T_vals,
      EB_res      = cnf_res,
      cnf_model   = cnf_res$cnf_model,
      adj_warning = 0
    )
  }else{

  # ---- KDE branch (default) ----
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
                                           verbose = verbose,
                                           sp_weights = sp_weights,
                                           monotone = monotone,
                                           monotone_step = monotone_step)
  lfdr_vec <- as.numeric(EB_res$lfdr)
  gamma <- 1-lfdr_vec
  adj_warning <- EB_res$adj_warning
  
  ret <- list(gamma = gamma,
              lfdr = lfdr_vec,
              pi0 = EB_res$pi0,
              p_vals = p_vals,
              T_vals = T_vals,
              EB_res = EB_res,
              adj_warning = adj_warning)
  }

  return(ret)
}





spatial_nonparam_kde_truncnorm <- function(T_vals,
                                           p_vals,
                                           Dfull_grid,
                                           neighbor_list, # list of neighbors
                                           h_spatial = 0.1,
                                           h_z = 1,          # KDE bandwidth in z multiplyer
                                           radius = NULL,       # neighbor radius in grid units (if NULL use 6)
                                           null_frac = 0.70,    # fraction of central Z to treat as local null
                                           lambda = 0.5,        # Storey tuning
                                           min_neighbors = 10,  # fallback if too few neighbors
                                           lambda_tau = 0.1,     # Penalty for estimate of sigma
                                           df = 40,
                                           H = NULL,            # True value of the signals
                                           eps = 1e-12,
                                           coord_leng,
                                           verbose = TRUE,
                                           sp_weights = NULL,   # output of build_sp_weights()
                                           monotone = FALSE,    # make each location's lfdr non-increasing in |T|
                                           monotone_step = 0.1) {

  n <- length(T_vals)
  if (is.null(h_z)) h_z <- 1
  # spatial weights do not change across iterations; sp_probe() builds them once
  if (is.null(sp_weights)) sp_weights <- build_sp_weights(neighbor_list, Dfull_grid, h_spatial, min_neighbors)
  w   <- sp_weights$w      # weight of each (location, neighbor) pair
  loc <- sp_weights$loc    # location index of each pair
  nb  <- sp_weights$nb     # neighbor index of each pair
  sw  <- sp_weights$sw     # sum of weights for each location
  # sum over each location's neighbors (pairs are sorted by location)
  loc_sum <- function(x) loc_sum_sorted(x, sp_weights$ends)

  # precompute p-values and global central threshold
  z_thres <- as.numeric(stats::quantile(abs(T_vals), probs = null_frac, na.rm = TRUE))
  z_abs_q <- abs(T_vals) <= z_thres

  global_sigma <- wt_truncnorm_sigma_ge1(T_vals[z_abs_q], rep(1,sum(z_abs_q)), t = z_thres)

  T_nb <- T_vals[nb]

  # ---- local bandwidth: weighted Silverman (see weighted_h) ----
  mu_loc  <- loc_sum(w * T_nb) / sw
  sd_loc  <- sqrt(loc_sum(w * (T_nb - mu_loc[loc])^2) / sw)
  h_z_vec <- 1.06 * pmax(sd_loc, 1) * sp_weights$n_nb^(-1/5)
  h_z_vec[!is.finite(h_z_vec) | h_z_vec <= 0] <- 0.3
  h_z_vec <- h_z * h_z_vec

  # ---- marginal f: weighted KDE on observed T in neighborhood ----
  f_hat_vec <- loc_sum(w * dnorm((T_vals[loc] - T_nb) / h_z_vec[loc])) / (h_z_vec * sw)
  f_hat_vec[!is.finite(f_hat_vec) | f_hat_vec <= 0] <- eps
  # KDE at 0, for the bandwidth warning
  f_at_0 <- loc_sum(w * dnorm(T_nb / h_z_vec[loc])) / (h_z_vec * sw)
  f_at_0[!is.finite(f_at_0) | f_at_0 <= 0] <- eps

  # ---- null f0 estimation ----
  # The weighted truncated-normal likelihood only depends on the data through
  # S0 = sum(w) and S2 = sum(w * T^2) over the central neighbors, so all local
  # sigmas are fit at once (see wt_truncnorm_sigma_ge1 for the single fit).
  cen_nb <- z_abs_q[nb]
  n_cen  <- loc_sum(1 * cen_nb)
  S0     <- loc_sum(w * cen_nb)
  S2     <- loc_sum(w * cen_nb * T_nb^2)
  # too few central points: use all neighbors
  all_nb <- n_cen < min_neighbors
  if (any(all_nb)) {
    S0[all_nb] <- sw[all_nb]
    S2[all_nb] <- loc_sum(w * T_nb^2)[all_nb]
  }
  n_null <- ifelse(all_nb, sp_weights$n_nb, n_cen)
  sigma_vec <- wt_truncnorm_sigma_ge1_vec(S0, S2, t = z_thres, lambda_tau = lambda_tau)
  # fallback: use global estimate of sigma
  sigma_vec[n_null < 3] <- global_sigma
  f0_hat_vec <- dnorm(T_vals, mean = 0, sd = sigma_vec)
  f0_hat_vec[!is.finite(f0_hat_vec) | f0_hat_vec <= 0] <- eps

  # ---- local pi0 via weighted Storey ----
  pi0_vec <- pmin(1, loc_sum(w * (p_vals[nb] > lambda)) / pmax((1 - lambda) * sw, eps))

  # ---- lfdr plug-in and clamp ----
  lfdr_vec <- pmax(0, pmin(1, pi0_vec * f0_hat_vec / f_hat_vec))
  lfdr_raw <- lfdr_vec
  if (monotone) {
    lfdr_vec <- monotone_lfdr_sp(T_vals, lfdr_vec, pi0_vec, sigma_vec, h_z_vec, sp_weights,
                                 step = monotone_step, eps = eps)
  }
  lfdr_at_0_vec <- f_at_0 / dnorm(0, mean = 0, sd = sigma_vec)

  ## Bandwidth warning:
  adj_warning = 0
  if(any(lfdr_at_0_vec > 1.5)){
    adj_warning = 1
  }

  EB_res <- list(lfdr = lfdr_vec,
                 lfdr_raw = lfdr_raw,   # before the monotone step (equal to lfdr if monotone = FALSE)
                 f = f_hat_vec,
                 f0 = f0_hat_vec,
                 pi0 = pi0_vec,
                 h_z = h_z_vec,
                 sigma = sigma_vec,
                 f_f0_rat_at_0_vec = lfdr_at_0_vec,
                 params = list(h_spatial = h_spatial, h_z = h_z, radius = radius, null_frac = null_frac,
                               monotone = monotone, monotone_step = monotone_step),
                 neighbor_list = sp_weights$neighbor_list,
                 adj_warning = adj_warning)

  if (verbose && length(coord_leng) == 2) {
    plot_e_step_sp(EB_res, T_vals, p_vals, Dfull_grid, coord_leng, H = H)
  }

  return(EB_res)
}


# ------------------------------------------------------------------
# build_sp_weights
#
# Spatial Gaussian weights exp(-d^2 / (2 h_spatial^2)) for every
# (location, neighbor) pair, stored as flat vectors sorted by location.
# Locations with fewer than min_neighbors neighbors use the full grid.
# Computed once per sp_probe() fit, since the weights do not change
# across iterations.
#
# Returns a list with w, loc, nb (one entry per pair), sw (sum of weights
# per location), n_nb (neighbors per location) and the neighbor_list
# actually used.
# ------------------------------------------------------------------
build_sp_weights <- function(neighbor_list, Dfull_grid, h_spatial, min_neighbors = 10) {
  n  <- length(neighbor_list)
  nl <- lapply(neighbor_list, function(v) if (length(v) < min_neighbors) seq_len(n) else v)
  n_nb <- lengths(nl)
  loc  <- rep.int(seq_len(n), n_nb)
  nb   <- unlist(nl, use.names = FALSE)
  w    <- exp(-Dfull_grid[cbind(nb, loc)]^2 / (2 * h_spatial^2))
  ends <- cumsum(n_nb)
  sw   <- loc_sum_sorted(w, ends)
  list(w = w, loc = loc, nb = nb, sw = sw, n_nb = n_nb, ends = ends, neighbor_list = nl)
}

# ------------------------------------------------------------------
# monotone_lfdr_sp
#
# Location-specific version of the monotone step in lfdr_t_func(). There,
# all locations share one lfdr curve, so sorting the observed T and taking
# cummax from the most extreme T inward gives
#   lfdr(T_j) <- max{ lfdr(t) : T_j <= t <= max(T) }   for T_j >= 0
# (and the mirror image for T_j < 0). Here each location j has its own
# curve lfdr_j(t) = min(1, pi0_j * phi(t / sigma_j) / sigma_j / f_j(t)), so
# the maximum is taken over location j's curve for |T_j| <= |t| <= t_max_j
# on the same side as T_j, where t_max_j is the largest |T| among j's
# neighbors. Capping at t_max_j keeps the tails of the local KDE, where
# there is no local data, from affecting the result. The curve is
# evaluated at grid points (spacing step) inside the window and at t_max_j
# itself. Locations whose running maximum has reached 1 are dropped, since
# they cannot change.
# ------------------------------------------------------------------
monotone_lfdr_sp <- function(T_vals, lfdr, pi0, sigma, h, sp_weights, step = 0.1, eps = 1e-12) {
  w    <- sp_weights$w
  loc  <- sp_weights$loc
  T_nb <- T_vals[sp_weights$nb]
  n_nb <- sp_weights$n_nb
  sw   <- sp_weights$sw
  a     <- abs(T_vals)
  t_max <- vapply(sp_weights$neighbor_list, function(v) max(a[v]), numeric(1))
  sgn   <- ifelse(T_vals >= 0, 1, -1)
  out   <- lfdr

  # lfdr_j at t = sgn_j * t_j for the active locations
  update <- function(act, t_act) {
    pr  <- which(act[loc])                    # pairs of the active locations
    lp  <- loc[pr]
    f_t <- loc_sum_sorted(w[pr] * dnorm((sgn[lp] * t_act[lp] - T_nb[pr]) / h[lp]),
                          cumsum(n_nb[act])) / (h[act] * sw[act])
    pmin(1, pi0[act] * dnorm(t_act[act], 0, sigma[act]) / pmax(f_t, eps))
  }

  for (g in seq(step, max(t_max), by = step)) {
    act <- a < g & g < t_max & out < 1
    if (!any(act)) next
    out[act] <- pmax(out[act], update(act, rep(g, length(a))))
  }
  act <- a < t_max & out < 1                  # right end of each window
  if (any(act)) out[act] <- pmax(out[act], update(act, t_max))
  out
}

# Group sums of x over consecutive blocks ending at positions ends (pairs
# sorted by location). Faster than rowsum(), which re-derives the groups.
loc_sum_sorted <- function(x, ends) {
  cs <- cumsum(x)[ends]
  cs - c(0, cs[-length(cs)])
}


# ------------------------------------------------------------------
# wt_truncnorm_sigma_ge1_vec
#
# Vectorized wt_truncnorm_sigma_ge1(): for each location, minimizes
#   S2 / (2 sigma^2) + S0 * (log(sigma) + log(Phi(t/sigma) - Phi(-t/sigma)))
#     + 0.5 * lambda_tau * (sigma - 1)^2
# over sigma in [lower, upper], where S0 = sum(w) and S2 = sum(w * z^2)
# over the local central values (the weighted truncated-normal negative
# log-likelihood, up to a constant). Uses golden-section search run in
# parallel over locations; n_iter = 50 gives an interval width of about
# (upper - lower) * 0.618^50 = 7e-10.
# ------------------------------------------------------------------
wt_truncnorm_sigma_ge1_vec <- function(S0, S2, t = 2, lambda_tau = 0.1,
                                       lower = 1, upper = 20, n_iter = 50) {
  nll <- function(s) {
    Znorm <- pnorm(t / s) - pnorm(-t / s)
    out <- S2 / (2 * s^2) + S0 * (log(s) + log(Znorm)) + 0.5 * lambda_tau * (s - 1)^2
    out[Znorm <= 1e-12] <- 1e20
    out
  }
  g  <- (sqrt(5) - 1) / 2
  a  <- rep(lower, length(S0)); b <- rep(upper, length(S0))
  x1 <- b - g * (b - a); x2 <- a + g * (b - a)
  f1 <- nll(x1); f2 <- nll(x2)
  for (it in seq_len(n_iter)) {
    # minimum in [a, x2] if f1 < f2, else in [x1, b]; one new point per location
    lft <- f1 < f2
    b[lft]  <- x2[lft];  a[!lft] <- x1[!lft]
    x2[lft] <- x1[lft];  f2[lft] <- f1[lft]
    x1[!lft] <- x2[!lft]; f1[!lft] <- f2[!lft]
    x_new <- ifelse(lft, b - g * (b - a), a + g * (b - a))
    f_new <- nll(x_new)
    x1[lft] <- x_new[lft];   f1[lft] <- f_new[lft]
    x2[!lft] <- x_new[!lft]; f2[!lft] <- f_new[!lft]
  }
  (a + b) / 2
}


# ------------------------------------------------------------------
# plot_e_step_sp
#
# 3 x 3 diagnostic panel for the spatial KDE E-step: true mask (or p-value
# histogram), observed T, lfdr, f0, f, pi0, KDE bandwidth, local null SD,
# and the marginal / null / lfdr curves at the location with max |T|.
# Works from stored E-step output, so it can be drawn after the fit.
#
# Arguments:
#   EB_res     : output of spatial_nonparam_kde_truncnorm() (E_step$EB_res)
#   T_vals     : test statistics (E_step$T_vals)
#   p_vals     : p-values (E_step$p_vals)
#   Dfull_grid : M x M distance matrix
#   coord_leng : c(nrow, ncol) of the 2D grid
#   H          : optional true signal mask (nrow x ncol)
# ------------------------------------------------------------------
plot_e_step_sp <- function(EB_res, T_vals, p_vals, Dfull_grid, coord_leng, H = NULL) {

  to_mat <- function(v) matrix(v, nrow = coord_leng[1], ncol = coord_leng[2], byrow = FALSE)
  img    <- function(m, main) fields::image.plot(t(m)[, nrow(m):1], main = main)

  # Curves at the location with the largest |T|, rebuilt from the stored
  # neighbor list, bandwidths, local null SD and local pi0.
  mx_i      <- which.max(abs(T_vals))
  neigh     <- EB_res$neighbor_list[[mx_i]]
  w_spatial <- exp(-Dfull_grid[neigh, mx_i]^2 / (2 * EB_res$params$h_spatial^2))
  x         <- seq(-6, 6, 0.1)
  f_est_x   <- vapply(x, function(x0) weighted_kde_eval(T_vals[neigh], w_spatial, x0,
                                                         EB_res$h_z[mx_i]), numeric(1))
  f0_est_x     <- dnorm(x, mean = 0, sd = EB_res$sigma[mx_i])
  f0_est_x_pi0 <- EB_res$pi0[mx_i] * f0_est_x
  lfdr_est_x   <- pmin(1, f0_est_x_pi0 / f_est_x)

  op <- par(no.readonly = TRUE)
  on.exit(tryCatch(par(op), error = function(e) NULL))
  par(mfrow = c(3,3), mar = c(2,2,2,2)+2)
  if (!is.null(H)){
    img(H, "True H (mask)")
  }else{
    hist(p_vals, main = "Observed p histogram", breaks = 40, col = gray(0.9))
  }
  img(to_mat(T_vals),       "Observed T")
  img(to_mat(EB_res$lfdr),  "Estimated lfdr")
  img(to_mat(EB_res$f0),    "Estimated f0")
  img(to_mat(EB_res$f),     "Estimated f")
  img(to_mat(EB_res$pi0),   "Estimated pi0")
  img(to_mat(EB_res$h_z),   "Bandwidth")
  img(to_mat(EB_res$sigma), "Estimated Null Sigma")
  big_vec <- c(f0_est_x, f_est_x, lfdr_est_x)
  big_vec[big_vec > 1] <- 1
  plot(x, f_est_x, type = "l", ylim = range(big_vec), main = "Distributions at max|T|")
  lines(x, f0_est_x, col = 2)
  lines(x, f0_est_x_pi0, col = 3)
  lines(x, lfdr_est_x, col = 4)
  abline(v = T_vals[mx_i])
  invisible(NULL)
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
  n_eff <- length(w)
  temp_h_z <- 1.06 * sd_z * n_eff^(-1/5)
  if (temp_h_z <= 0 || is.na(temp_h_z)) temp_h_z <- 0.3
  temp_h_z
}

# Calculates a weighted standard deviation.
weighted_sd <- function(x, w) {
  mu <- sum(w * x) / sum(w)
  sqrt(sum(w * (x - mu)^2) / sum(w))
}


# ------------------------------------------------------------------
# spatial_lfdr_cnf
#
# Spatial lfdr estimation using a conditional normalizing flow (CNF)
# neural-network mixture-of-logistics model for the marginal density
# f(z | s), a local weighted Storey estimator for pi0(s), and an
# optional local null variance estimation via a truncated N(0, sigma^2)
# model with shrinkage.
#
# Requires the 'torch' package and the helper functions:
#   train_cnf(), estimate_f_cnf(), lfdr_threshold()
#
# Arguments
# ---------
# Z             : nrow x ncol matrix of observed Z-/T-statistics
# coords        : n x 2 coordinate matrix (column-major order, n = nrow*ncol)
# h_spatial     : spatial Gaussian kernel bandwidth; auto-computed from
#                 Dfull_grid if NULL
# Dfull_grid    : n x n distance matrix; used only when h_spatial = NULL
# radius        : neighbour radius in grid units (default: max(3, min_dim/8))
# null_frac     : central fraction of |Z| used to identify null obs when
#                 estimating sigma locally (default 0.70)
# lambda        : Storey tuning parameter for pi0 (default 0.5)
# min_neighbors : fall back to all locations if neighbourhood < this (10)
# alpha         : target FDR level for lfdr thresholding (default 0.1)
# estimate_sigma: if TRUE (default), estimate null sigma locally via
#                 weighted truncated-normal MLE with shrinkage.
#                 if FALSE, fix null sigma = sigma_null at every location.
# sigma_null    : fixed null SD used when estimate_sigma = FALSE (default 1)
# lambda_tau    : shrinkage penalty toward sigma = 1 in local MLE;
#                 only used when estimate_sigma = TRUE (default 0.001)
# H             : true signal indicator matrix (optional; for diagnostics)
# eps           : numerical floor (default 1e-12)
# cnf_model     : a previously fitted CNF model returned by an earlier call
#                 (the $cnf_model element). When supplied, training is skipped
#                 and all cnf_* tuning parameters below are ignored.
# cnf_K         : MoL mixture components (default 10)
# cnf_hidden_dim: conditioner network width (default 64)
# cnf_n_hidden  : conditioner hidden layers (default 2)
# cnf_n_epochs  : Adam training epochs (default 1000)
# cnf_lr        : Adam learning rate (default 1e-3)
# cnf_batch_size: mini-batch size for training (default 256)
# verbose       : print progress and plots (default TRUE)
#
# Returns a list with:
#   lfdr, f, f0, pi0, sigma_mat  — matrices (nrow x ncol)
#   thresh_info                  — output of lfdr_threshold()
#   cnf_model                    — fitted CNF model object
#   params                       — list of run parameters
#   neighbor_list                — list of neighbour index vectors
# ------------------------------------------------------------------
spatial_lfdr_cnf <- function(Z,
                             coords,
                             h_spatial      = NULL,
                             Dfull_grid     = NULL,
                             radius         = NULL,
                             null_frac      = 0.70,
                             lambda         = 0.5,
                             min_neighbors  = 10,
                             alpha          = 0.1,
                             estimate_sigma = TRUE,
                             sigma_null     = 1,
                             lambda_tau     = 0.001,
                             H              = NULL,
                             eps            = 1e-12,
                             cnf_model      = NULL,
                             cnf_K          = 10,
                             cnf_hidden_dim = 64,
                             cnf_n_hidden   = 2,
                             cnf_n_epochs   = 1000,
                             cnf_lr         = 1e-3,
                             cnf_batch_size = 256,
                             verbose        = TRUE) {

  nrow_ <- nrow(Z)
  ncol_ <- ncol(Z)
  n     <- nrow_ * ncol_
  Zvec  <- as.numeric(Z)

  # ---- radius and spatial bandwidth ----
  if (is.null(radius)) radius <- max(3, min(nrow_, ncol_) / 8)
  if (is.null(h_spatial)) {
    if (is.null(Dfull_grid))
      stop("Provide either h_spatial or Dfull_grid.")
    h_spatial <- median(Dfull_grid[upper.tri(Dfull_grid)]) / 6
    if (is.na(h_spatial) || h_spatial <= 0) h_spatial <- 3
  }

  if (verbose)
    message("spatial_lfdr_cnf: n = ", n,
            ", radius = ", radius,
            ", h_spatial = ", round(h_spatial, 3),
            ", estimate_sigma = ", estimate_sigma,
            if (!estimate_sigma) paste0("  (fixed sigma_null = ", sigma_null, ")") else "")

  # ---- neighbour lists ----
  neighbor_list <- vector("list", n)
  radius2 <- radius^2
  for (i in seq_len(n)) {
    d2 <- (coords[, 1] - coords[i, 1])^2 + (coords[, 2] - coords[i, 2])^2
    neighbor_list[[i]] <- which(d2 <= radius2)
  }

  # ---- global quantities ----
  pvals   <- 2 * (1 - pnorm(abs(Zvec)))
  z_thres <- as.numeric(quantile(abs(Zvec), probs = null_frac, na.rm = TRUE))
  z_abs_q <- abs(Zvec) <= z_thres   # central-null mask (logical, length n)

  global_sigma <- if (estimate_sigma) {
    wt_truncnorm_sigma_ge1(Zvec[z_abs_q], rep(1, sum(z_abs_q)), t = z_thres,
                           lambda_tau = lambda_tau)
  } else {
    sigma_null
  }

  if (verbose && estimate_sigma)
    message("  global null sigma = ", round(global_sigma, 3))

  # ---- train CNF (or reuse supplied model) and evaluate f ----
  if (is.null(cnf_model)) {
    if (verbose) message("Training CNF (neural network MoL, requires torch) ...")
    cnf_model <- train_cnf(Z, coords,
                           K          = cnf_K,
                           hidden_dim = cnf_hidden_dim,
                           n_hidden   = cnf_n_hidden,
                           n_epochs   = cnf_n_epochs,
                           lr         = cnf_lr,
                           batch_size = cnf_batch_size,
                           verbose    = verbose)
  } else {
    if (verbose) message("Using supplied CNF model (skipping training).")
  }
  f_hat_vec <- estimate_f_cnf(Z, coords, cnf_model)

  # ---- per-location loop: f0, pi0, lfdr ----
  f0_hat_vec <- numeric(n)
  pi0_vec    <- numeric(n)
  lfdr_vec   <- numeric(n)
  sigma_vec  <- numeric(n)

  for (i in seq_len(n)) {

    neigh <- neighbor_list[[i]]
    if (length(neigh) < min_neighbors) neigh <- seq_len(n)

    d2        <- (coords[neigh, 1] - coords[i, 1])^2 +
                 (coords[neigh, 2] - coords[i, 2])^2
    w_spatial <- exp(-d2 / (2 * h_spatial^2))
    if (sum(w_spatial) <= 0) w_spatial <- rep(1, length(w_spatial))

    # ---- null sigma ----
    if (estimate_sigma) {
      central_mask <- z_abs_q[neigh]
      if (sum(central_mask) < min_neighbors)
        central_mask <- rep(TRUE, length(neigh))
      samples_null  <- Zvec[neigh][central_mask]
      weights_null  <- w_spatial[central_mask]
      if (length(samples_null) < 3) {
        est_sigma <- global_sigma
      } else {
        est_sigma <- wt_truncnorm_sigma_ge1(samples_null, weights_null,
                                            t = z_thres,
                                            lambda_tau = lambda_tau)
      }
    } else {
      est_sigma <- sigma_null
    }
    sigma_vec[i] <- est_sigma

    # ---- f0 ----
    f0_i <- dnorm(Zvec[i], 0, est_sigma)
    if (!is.finite(f0_i) || f0_i <= 0) f0_i <- eps
    f0_hat_vec[i] <- f0_i

    # ---- pi0 via local Storey ----
    numer      <- sum(w_spatial * (pvals[neigh] > lambda))
    denom      <- (1 - lambda) * sum(w_spatial)
    pi0_vec[i] <- min(1, numer / max(denom, eps))

    # ---- lfdr ----
    lfdr_vec[i] <- max(0, min(1, pi0_vec[i] * f0_i / max(f_hat_vec[i], eps)))
  }

  # ---- reshape to matrices ----
  lfdr_mat  <- matrix(lfdr_vec,   nrow = nrow_, ncol = ncol_, byrow = FALSE)
  f_mat     <- matrix(f_hat_vec,  nrow = nrow_, ncol = ncol_, byrow = FALSE)
  f0_mat    <- matrix(f0_hat_vec, nrow = nrow_, ncol = ncol_, byrow = FALSE)
  pi0_mat   <- matrix(pi0_vec,    nrow = nrow_, ncol = ncol_, byrow = FALSE)
  sigma_mat <- matrix(sigma_vec,  nrow = nrow_, ncol = ncol_, byrow = FALSE)

  # ---- diagnostic plots ----
  if (verbose) {
    op <- par(no.readonly = TRUE)
    par(mfrow = c(3, 3), mar = c(2, 2, 2, 1))
    if (!is.null(H)) {
      fields::image.plot(t(H)[, nrow(H):1],       main = "True H (mask)")
    } else {
      hist(pvals, main = "Observed p histogram", breaks = 40, col = gray(0.9))
    }
    fields::image.plot(t(Z)[, nrow(Z):1],          main = "Observed Z")
    fields::image.plot(t(lfdr_mat)[, nrow_:1],     main = "Estimated lfdr (CNF)")
    fields::image.plot(t(f0_mat)[, nrow_:1],       main = "Estimated f0")
    fields::image.plot(t(f_mat)[, nrow_:1],        main = "Estimated f (CNF)")
    fields::image.plot(t(pi0_mat)[, nrow_:1],      main = "Estimated pi0")
    if (estimate_sigma) {
      fields::image.plot(t(sigma_mat)[, nrow_:1],  main = "Estimated null sigma")
    } else {
      plot.new()
      title(main = paste0("Fixed null sigma = ", sigma_null))
    }
    plot.new()
    plot(sort(as.numeric(lfdr_mat)), type = "l", main = "Sorted lfdr values")
    tryCatch(par(op), error = function(e) NULL)
  }

  # ---- threshold ----
  res_thresh <- lfdr_threshold(lfdr_mat = lfdr_mat, alpha = alpha,
                               true_H = H, nrow_ = nrow_, ncol_ = ncol_,
                               verbose = verbose)

  if (verbose) {
    op2 <- par(no.readonly = TRUE)
    par(mfrow = c(1, 3))
    if (!is.null(H)) fields::image.plot(t(H)[, nrow_:1], main = "True Signals")
    fields::image.plot(t(lfdr_mat)[, nrow_:1],   main = "Estimated lfdr")
    fields::image.plot(t(res_thresh$reject_mat)[, nrow_:1],
                       main = paste0("Rejections (alpha = ", alpha, ")"))
    tryCatch(par(op2), error = function(e) NULL)
  }

  return(list(
    lfdr         = lfdr_mat,
    f            = f_mat,
    f0           = f0_mat,
    pi0          = pi0_mat,
    sigma_mat    = sigma_mat,
    thresh_info  = res_thresh,
    cnf_model    = cnf_model,
    params       = list(
      h_spatial      = h_spatial,
      radius         = radius,
      null_frac      = null_frac,
      estimate_sigma = estimate_sigma,
      sigma_null     = sigma_null,
      lambda_tau     = lambda_tau,
      cnf_K          = cnf_K,
      cnf_hidden_dim = cnf_hidden_dim,
      cnf_n_hidden   = cnf_n_hidden,
      cnf_n_epochs   = cnf_n_epochs,
      cnf_lr         = cnf_lr
    ),
    neighbor_list = neighbor_list
  ))
}
