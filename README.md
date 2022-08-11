
# Introduction

This repository contains the R software tools to run a UNivariate
approach to HIgh-Dimensional linear regression via the EM algorithm
(UNHIDEM). UNHIDEM is an innovative approach to performing
high-dimensional linear regression that provide covariate-wide
estimation through a focus on the adjustment approximation. Our approach
considers the unknown adjustment predictor as missing data and estimates
the regression parameter using a expectation-maximization (EM) approach.

The framework used by UNHIDEM has similarities to high-dimensional
causal inference in the focus (for each predictor) on estimating the
unknown adjustment variable and to Bayesian variable selection through
the use of latent binary variables indicating a predictor is non-null
(Ročková and George 2014). The estimation approach is similar to an
application of the partitioned Expectation–Conditional Maximization
(PECM, Meng and Rubin 1992, 1993) algorithm on the complete data
likelihood, where for each predictor the parameters have been expanded
similar to the parameter-expanded EM (PX-EM, Liu, Rubin, and Wu 1998).

# Example of UNHIDEM

The following is demonstration of how to implement UNHIDEM with
simulated data. Please note that the package is still a work in progress
and we are continually making improvements.

Install the package.

``` r
library(devtools)
install_github(repo="alexmclain/UNHIDEM", subdir="unhidem")
```

Load the package and the data.

``` r
library(unhidem)
data(Sim_data)
attach(Sim_data) 
M <- dim(Z)[2] 
M1 <- length(signal) 
```

Here, **Sim\_data** contains the following elements:

-   **Y**: vector of outcome variables for the training set,
-   **Z**: *n* × *M* covariate matrix for the training set,
-   **beta\_tr**: true value of beta coefficients,
-   **signal**: indicies of the *M*<sub>1</sub> non-null predictors,
-   **Y\_test**: outcome variable for the test set, and
-   **Z\_test**: covariate matrix for the test set

where *n* = 400, *M* = 10, 000 and *M*<sub>1</sub> = 100.

Set convergence criteria, alpha, and plot indicator:

``` r
ep <- 0.004
alpha <- 0.05
plot_ind <- TRUE
```

In this first run of the analysis we include the true signal
(*η*<sub>*i*</sub> = *Z*<sub>*i*</sub>*β*) and indices of the non-null
beta coefficients (signal). This information will be used to create a
plot of the MSE and the number of rejections with an indication of
whether FDR ≤ *α* by iteration (*t*).

``` r
eta_i <- apply(t(Z)*beta_tr,2,sum) 
full_res <- unhidem( Y = Y, Z = Z, alpha = alpha, ep = ep, plot_ind = 
                       plot_ind, eta_i = eta_i, signal = signal)
```

An example of performing prediction (of training data).

``` r
pred_res <- predict_unhidem_func( full_res, Z = Z, alpha = alpha)
head(pred_res)
```

Applying multiple testing corrections to results

``` r
E_step_final <- full_res$E_step
MT_res <- mtr_func(E_step_final, alpha, signal)
```

Summary of multiple testing methods:

``` r
MT_res$Bonf_sum
MT_res$BH_sum
## Confusion matrix for BH procedure
table((1:M %in% signal),MT_res$BH_res$BH)
```

(note: the ability of UNHIDEM to appropriately control the FDR has not
been studied).

In this second run we’ll include test data instead of *η*<sub>*i*</sub>
and *signal*. This run will give identical results to the first, but the
plot will show the test MSE and the total number of rejections by
iteration.

``` r
dim(Z_test)
length(Y_test)

full_res <- unhidem( Y = Y, Z = Z, alpha = alpha, ep = ep, plot_ind = 
                       plot_ind, Y_test = Y_test, Z_test = Z_test)
```

Predict for the test observations

``` r
pred_res_test <- predict_unhidem_func(full_res, Z = Z_test, 
                                       alpha = alpha)

## The true signal for test data
eta_test <- apply(t(Z_test)*beta_tr,2,sum) 
```

Proportion of test CIs and PIs that contain the true signal and the
(unobserved) test observation, respectively.

``` r
ecp_CI <- mean(1*I(eta_test>pred_res_test$CI_L & eta_test<pred_res_test$CI_U))
ecp_PI <- mean(1*I(Y_test>pred_res_test$PI_L & Y_test<pred_res_test$PI_U))
c(ecp_CI, ecp_PI)
detach(Sim_data)
```

# Example of UNHIDEM with covariate data

This example will again fit a high-dimensional linear regression model,
but this time additional covariates that are not subjected to the
“sparsity” assumption are included. That is, we fit the model

*Y*<sub>*i*</sub> = *μ* + *Z*<sub>*i*</sub>*β* + *X*<sub>*i*</sub>*β*<sub>*X*</sub> + *ϵ*

where *Z* has dimension *M* ≫ *n* and *X* has dimension *p* &lt; *n*
(here *p* = 2).

Load in the data:

``` r
data(Sim_data_cov)
attach(Sim_data_cov) 
M <- dim(Z)[2] 
M1 <- length(signal) 
eta_i <- apply(t(Z)*beta_tr,2,sum) 
```

Here, **Sim\_data** contains the following elements:

-   **Y**: vector of outcome variables for the training set,
-   **Z**: *n* × *M* predictor matrix for the training set,
-   **X**: *n* × 2 covariate matrix for the training set (not subjected
    to the “sparsity” assumption),
-   **beta\_tr**: true value of *β*,
-   **beta\_X\_tr**: true value of *β*<sub>*X*</sub>,
-   **signal**: indicies of the *M*<sub>1</sub> non-null predictors,

where *n* = 400, *M* = 10000 and *M*<sub>1</sub> = 100. The true signal
term, *η*<sub>*i*</sub>, only includes the impact of **Z** (not of
**X**).

In this analysis we include *η*<sub>*i*</sub> and **signal** for
plotting purposes.

``` r
full_res <- unhidem(Y = Y, Z = Z, X = X, alpha = alpha, ep = ep, plot_ind 
                    = plot_ind, eta_i = eta_i, signal = signal)
```

Total number of iterations

``` r
full_res$count
```

Final estimates of the impact of X versus the true values:

``` r
data.frame(true_values = beta_X_tr, full_res$Calb_mod$res_data[-2,])
```

Compare to a standard linear model of *X* on *Y*:

``` r
summary(lm(Y~X$Cont_cov + X$Binary_cov))$coefficients
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Liuetal98" class="csl-entry">

Liu, Chuanhai, Donald B. Rubin, and Ying Nian Wu. 1998. “<span
class="nocase">Parameter expansion to accelerate EM: The PX-EM
algorithm</span>.” *Biometrika* 85 (4): 755–70.
<https://doi.org/10.1093/biomet/85.4.755>.

</div>

<div id="ref-MenRub92" class="csl-entry">

Meng, Xiao-Li, and Donald B. Rubin. 1992. “Recent Extensions to the EM
Algorithm.” In *Bayesian Statistics, 4 (Peñı́scola, 1991)*, 307–20.
Oxford Univ. Press, New York.

</div>

<div id="ref-MenRub93" class="csl-entry">

———. 1993. “<span class="nocase">Maximum likelihood estimation via the
ECM algorithm: A general framework</span>.” *Biometrika* 80 (2): 267–78.
<https://doi.org/10.1093/biomet/80.2.267>.

</div>

<div id="ref-VerGeo14" class="csl-entry">

Ročková, Veronika, and Edward I. George. 2014. “EMVS: The EM Approach to
Bayesian Variable Selection.” *Journal of the American Statistical
Association* 109 (506): 828–46.
<https://doi.org/10.1080/01621459.2013.869223>.

</div>

</div>
