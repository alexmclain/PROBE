generate_factor_data <- function(n, p, s, SNR, K = 10, n_test = 1000) {
  # --- Design matrix via latent factor model ---
  Z <- matrix(rnorm(n * K), n, K)
  Lambda <- matrix(rnorm(p * K), p, K)
  E <- matrix(rnorm(n * p), n, p)
  X_raw <- Z %*% t(Lambda) + E
  
  # --- Center and scale X column-wise ---
  X <- scale(X_raw, center = TRUE, scale = TRUE)
  col_means <- attr(X, "scaled:center")
  col_sds <- attr(X, "scaled:scale")
  X <- as.matrix(X)
  
  # --- Sparse coefficient vector ---
  signal_idx <- sort(sample(seq_len(p), s))
  beta_true <- rep(0, p)
  beta_true[signal_idx] <- runif(s, 0, 4)
  
  # --- Response, calibrated to target SNR ---
  eta_true <- X %*% beta_true
  sigma_eps <- sqrt(var(as.numeric(eta_true)) / SNR)
  Y <- as.numeric(eta_true + rnorm(n, sd = sigma_eps))
  Y <- Y - mean(Y)
  
  # --- Independent test set (same scaling) ---
  Z_test <- matrix(rnorm(n_test * K), n_test, K)
  E_test <- matrix(rnorm(n_test * p), n_test, p)
  X_test_raw <- Z_test %*% t(Lambda) + E_test
  X_test <- as.matrix(scale(X_test_raw, center = col_means, scale = col_sds))
  Y_test <- as.numeric(X_test %*% beta_true + rnorm(n_test, sd = sigma_eps))
  
  list(X = X, Y = Y, X_test = X_test, Y_test = Y_test,
       beta_true = beta_true, sigma_eps = sigma_eps)
}