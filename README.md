
# Introduction

This repository contains the R software tools to run the PaRtitiOned
empirical Bayes Ecm (PROBE) algorithm. We give a brief explanation of
the method below. For more details see McLain, Zgodic, and Bondell
(2022).

Minimal prior assumptions on the parameters are used through the use of
plug-in empirical Bayes estimates of hyperparameters. Efficient maximum
(MAP) estimation is completed through a Parameter-Expanded
Expectation-Conditional-Maximization (PX-ECM) algorithm. The PX-ECM
results in a robust computationally efficient coordinate-wise
optimization, which adjusts for the impact of other predictor variables.
The completion of the E-step uses an approach motivated by the popular
two-groups approach to multiple testing. The PROBE algorithm is applied
to sparse high-dimensional linear regression, which can be completed
using one-at-a-time or all-at-once type optimization. PROBE is a novel
alternative to Markov chain Monte Carlo (Liang et al. 2008; Bondell and
Reich 2012 ; Chae, Lin, and Dunson 2019), empirical Bayes (George and
Foster 2000; Martin, Mess, and Walker 2017; Martin and Tang 2020), and
Variational Bayes (Carbonetto and Stephens 2012; Blei, Kucukelbir, and
McAuliffe 2017; Ray and Szabó 2021) approaches to fitting sparse linear
models.

Our proposed method performs Bayesian variable selection with an
uninformative spike-and-slab prior on the regression parameters, which
has not yet been used in the high-dimensional setting. We focus on MAP
estimation of regression parameters, latent variable selection
indicators, and the residual variance. We use a quasi Parameter-Expanded
Expectation-Conditional-Maximization (PX-ECM) which is a combination of
the ECM Dyk, Meng, and Rubin (1995) and PX-EM (Liu, Rubin, and Wu 1998).
With the standard EM, the M-step would require optimization with respect
to the high-dimensional regression parameter, which is not feasible with
uninformative priors. Here, the ECM is used to break the M-step into a
coordinate-wise optimization procedure. The PX portion of the algorithm
adds a parameter which scales the impact of the remaining predictor
variables. Unlike other coordinate-wise optimization approaches, this
does not assume the impact of other predictor variables is fixed. The
benefits of this formulation are that the impact of the remaining
predictor variables only needs to be known up to a multiplicative
constant, we can account for the increase in variability due to the
dependence between the predictor and the impact of the remaining
predictor variables, and it adds stability to the algorithm. The E-step
consists of estimating the probability of the latent variable selection
indicators via updating hyperparameters with plug-in empirical Bayes
estimates, which is motivated by the popular two-groups approach to
multiple testing (Efron et al. 2001; Sun and Cai 2007).

One-at-a-time and all-at-once variants of the PROBE algorithm have been
developed. The examples below focus on the all-at-once variant.

# Example of PROBE

The following is demonstration of how to implement all-at-once PROBE
with simulated data. Please note that the package is still a work in
progress and we are continually making improvements.

Install the package.

``` r
library(devtools)
install_github(repo="alexmclain/PROBE", subdir="probe")
```

Load the package and the data.

``` r
library(probe)
data(Sim_data)
attach(Sim_data) 
M <- dim(X)[2] 
M1 <- length(signal) 
```

Here, **Sim_data** contains the following elements:

- **Y**: vector of outcome variables for the training set,
- **X**: $n \times M$ covariate matrix for the training set,
- **beta_tr**: true value of beta coefficients, and
- **signal**: indicies of the $M_1$ non-null predictors

where $n=400$, $M=10000$ and $M_1=100$.

Set convergence criteria, alpha, and plot indicator:

``` r
ep <- 0.01
alpha <- 0.05
plot_ind <- TRUE
```

In this first run of the analysis we include the true signal
($\eta_i= X_i \beta$) and indices of the non-null beta coefficients
(signal). This information will be used to create a plot of the MSE and
the number of rejections with an indication of whether FDR$\leq \alpha$
by iteration ($t$).

``` r
eta_i <- apply(t(X)*beta_tr,2,sum) 
full_res <- probe( Y = Y, X = X, alpha = alpha, ep = ep, plot_ind = 
                       plot_ind, eta_i = eta_i, signal = signal)
```

An example of performing prediction (of training data).

``` r
pred_res <- predict_probe_func( full_res, X = X, alpha = alpha)
head(pred_res)
```

How the credible and prediction intervals are constructed is discussed
in (Zgodic et al. 2023).

### Re-run with test data

In this second run we’ll include test data instead of $\eta_i$ and
*signal*. This run will give identical results to the first, but the
plot will show the test MSE and the total number of rejections by
iteration.

Load in some test data.

``` r
data(Sim_data_test)
attach(Sim_data_test) 
dim(X_test)
length(Y_test)
```

Here, **Sim_data_test** contains the following elements:

- **Y_test**: outcome variable for the test set, and
- **X_test**: covariate matrix for the test set.

Run using the test data. This information will be used to create a plot
of the test MSE by iteration. This is for demonstration only. The test
data does not select the number of iterations (the algorithm stops once
the convergence criteria is met).

``` r
full_res <- probe( Y = Y, X = X, alpha = alpha, ep = ep, plot_ind = 
                       plot_ind, Y_test = Y_test, X_test = X_test)
```

Predict for the test observations

``` r
pred_res_test <- predict_probe_func(full_res, X = X_test, 
                                       alpha = alpha)

## The true signal for test data
eta_test <- apply(t(X_test)*beta_tr,2,sum) 
plot(pred_res_test$Pred, eta_test)
```

Proportion of test CIs and PIs that contain the true signal and the
(unobserved) test observation, respectively.

``` r
ecp_CI <- mean(1*I(eta_test>pred_res_test$CI_L & eta_test<pred_res_test$CI_U))
ecp_PI <- mean(1*I(Y_test>pred_res_test$PI_L & Y_test<pred_res_test$PI_U))
c(ecp_CI, ecp_PI)
detach(Sim_data)
```

See Zgodic et al. (2023) for a thorough examination of the coverage
probabilities of the prediction intervals.

# Example of PROBE with covariate data

This example will again fit a high-dimensional linear regression model,
but this time additional covariates that are not subjected to the
“sparsity” assumption are included. That is, we fit the model

$$Y_i = \mu +  X_i \beta + Z_i \beta_Z + \epsilon$$

where $X$ has dimension $M \gg n$ and $Z$ has dimension $p<n$ (here
$p=2$). How non-sparse covariates are brought into the model is
discussed in Zgodic et al. (2023).

Load in the data:

``` r
data(Sim_data_cov)
attach(Sim_data_cov) 
M <- dim(X)[2] 
M1 <- length(signal) 
eta_i <- apply(t(X)*beta_tr,2,sum) 
```

Here, **Sim_data** contains the following elements:

- **Y**: vector of outcome variables for the training set,
- **X**: $n \times M$ predictor matrix for the training set,
- **Z**: $n \times 2$ covariate matrix for the training set (not
  subjected to the “sparsity” assumption),
- **beta_tr**: true value of $\beta$,
- **beta_Z_tr**: true value of $\beta_Z$,
- **signal**: indicies of the $M_1$ non-null predictors,

where $n=400$, $M=10000$ and $M_1=100$. The true signal term, $\eta_i$,
only includes the impact of **X** (not of **Z**).

In this analysis we include $\eta_i$ and **signal** for plotting
purposes.

``` r
full_res <- probe(Y = Y, X = X, Z = Z, alpha = alpha, ep = ep, plot_ind 
                    = plot_ind, eta_i = eta_i, signal = signal)
```

Total number of iterations

``` r
full_res$count
```

Final estimates of the impact of Z versus the true values:

``` r
data.frame(true_values = beta_Z_tr, full_res$Calb_mod$res_data[-2,])
```

Compare to a standard linear model of $Z$ on $Y$:

``` r
summary(lm(Y~Z$Cont_cov + Z$Binary_cov))$coefficients
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Bleetal17" class="csl-entry">

Blei, David M, Alp Kucukelbir, and Jon D McAuliffe. 2017. “Variational
Inference: A Review for Statisticians.” *Journal of the American
Statistical Association* 112 (518): 859–77.

</div>

<div id="ref-BonRei12" class="csl-entry">

Bondell, Howard D, and Brian J Reich. 2012. “Consistent High-Dimensional
Bayesian Variable Selection via Penalized Credible Regions.” *Journal of
the American Statistical Association* 107 (500): 1610–24.

</div>

<div id="ref-CarSte12" class="csl-entry">

Carbonetto, Peter, and Matthew Stephens. 2012. “Scalable Variational
Inference for Bayesian Variable Selection in Regression, and Its
Accuracy in Genetic Association Studies.” *Bayesian Analysis* 7 (1):
73–108. <https://doi.org/10.1214/12-BA703>.

</div>

<div id="ref-Chaeetal19" class="csl-entry">

Chae, Minwoo, Lizhen Lin, and David B Dunson. 2019. “Bayesian Sparse
Linear Regression with Unknown Symmetric Error.” *Information and
Inference: A Journal of the IMA* 8 (3): 621–53.

</div>

<div id="ref-vanetal95" class="csl-entry">

Dyk, David A. van, Xiao-Li Meng, and Donald B. Rubin. 1995. “Maximum
Likelihood Estimation via the ECM Algorithm: Computing the Asymptotic
Variance.” *Statist. Sinica* 5 (1): 55–75.

</div>

<div id="ref-Efron2001" class="csl-entry">

Efron, Bradley, Robert Tibshirani, John D. Storey, and Virginia Tusher.
2001. “Empirical Bayes Analysis of a Microarray Experiment.” *J. Amer.
Statist. Assoc.* 96 (456): 1151–60.
<https://doi.org/10.1198/016214501753382129>.

</div>

<div id="ref-GeoFos00" class="csl-entry">

George, Edward I, and Dean P Foster. 2000. “Calibration and Empirical
Bayes Variable Selection.” *Biometrika* 87 (4): 731–47.

</div>

<div id="ref-Lietal08" class="csl-entry">

Liang, Feng, Rui Paulo, German Molina, Merlise A Clyde, and Jim O
Berger. 2008. “Mixtures of g Priors for Bayesian Variable Selection.”
*Journal of the American Statistical Association* 103 (481): 410–23.

</div>

<div id="ref-Liuetal98" class="csl-entry">

Liu, Chuanhai, Donald B. Rubin, and Ying Nian Wu. 1998.
“<span class="nocase">Parameter expansion to accelerate EM: The PX-EM
algorithm</span>.” *Biometrika* 85 (4): 755–70.
<https://doi.org/10.1093/biomet/85.4.755>.

</div>

<div id="ref-Maretal17" class="csl-entry">

Martin, Ryan, Raymond Mess, and Stephen G Walker. 2017. “Empirical Bayes
Posterior Concentration in Sparse High-Dimensional Linear Models.”
*Bernoulli* 23 (3): 1822–47.

</div>

<div id="ref-MarTan20" class="csl-entry">

Martin, Ryan, and Yiqi Tang. 2020. “Empirical Priors for Prediction in
Sparse High-Dimensional Linear Regression.” *J. Mach. Learn. Res.* 21:
1–30.

</div>

<div id="ref-McLetal22" class="csl-entry">

McLain, Alexander C., Anja Zgodic, and Howard Bondell. 2022. “Sparse
High-Dimensional Linear Regression with a Partitioned Empirical Bayes
ECM Algorithm.” <https://arxiv.org/abs/2209.08139>.

</div>

<div id="ref-MenRub92" class="csl-entry">

Meng, Xiao-Li, and Donald B. Rubin. 1992. “Recent Extensions to the EM
Algorithm (with Discussion).” In *Bayesian Statistics 4*, edited by J.
M. Bernardo, J. O. Berger, A. P. Dawid, and A. F. M. Smith, 307–20.
Oxford University Press.

</div>

<div id="ref-MenRub93" class="csl-entry">

———. 1993. “<span class="nocase">Maximum likelihood estimation via the
ECM algorithm: A general framework</span>.” *Biometrika* 80 (2): 267–78.
<https://doi.org/10.1093/biomet/80.2.267>.

</div>

<div id="ref-KolBot21" class="csl-entry">

Ray, Kolyan, and Botond Szabó. 2021. “Variational Bayes for
High-Dimensional Linear Regression with Sparse Priors.” *Journal of the
American Statistical Association* 0 (0): 1–12.
<https://doi.org/10.1080/01621459.2020.1847121>.

</div>

<div id="ref-Sun2007" class="csl-entry">

Sun, Wenguang, and T. Tony Cai. 2007. “Oracle and Adaptive Compound
Decision Rules for False Discovery Rate Control.” *J. Amer. Statist.
Assoc.* 102 (479): 901–12. <https://doi.org/10.1198/016214507000000545>.

</div>

<div id="ref-Zgoetal23" class="csl-entry">

Zgodic, Anja, Ray Bai, Jiajia Zhang, Yuan Wang, Chris Rorden, and
Alexander McLain. 2023. “Heteroscedastic Sparse High-Dimensional Linear
Regression with a Partitioned Empirical Bayes ECM Algorithm.”
<https://arxiv.org/abs/2309.08783>.

</div>

</div>
