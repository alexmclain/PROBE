
# Introduction

The code in this folder can be use to reproduce the simulations given in
the paper. The outcome was generated via
$$Y_i = \mathbf{X}_i\mathbf{\gamma}\mathbf{\beta} + \mathbf{\epsilon}_i $$
where $\mathbf{\beta} \sim \mathcal{U}(0,2\eta)$, and
$\epsilon_i \sim N(0,\sigma^2)$. To generate $\mathbf{X}$ and
$\mathbf{\gamma}$, we considered (i) squared exponential and (ii)
grouped block diagonal covariance structures.

For the squared exponential covariance, each $\mathbf{X}_m$ has a two
dimensional coordinate
$\mathbf{d}_m = (d_{1m},d_{2m}) \in (1,\ldots,\sqrt{M})^2$; continuous
and binary $\mathbf{X}$ were considered. For continuous $\mathbf{X}$, we
used $\mathbf{X}_i \sim^{iid} MVN(0,\Sigma_{s})$ where
$\Sigma_{s}(\mathbf{d}_m,\mathbf{d}_{m^\prime}) = \exp\{-\lVert (\mathbf{d}_m - \mathbf{d}_{m^\prime})/s \rVert_2^2\}$
and $s=10$. The binary setting used
$\mathbf{X}^b_{i} = I( \mathbf{X}_{i}< 0)$ where $I(\cdot)$ denotes the
indicator function. For $\mathbf{\gamma}$, we set
$\gamma_m = I\{A_m<Q_{\pi}\}$ where
$(A_1,\ldots,A_M) \sim MVN(0,\Sigma_{20})$ and $Q_{\pi}$ was such that
$M_1/M=\pi$. Section C of the Supplementary Materials presents examples
of $\mathbf{\gamma}\mathbf{\beta}$ and $\mathbf{X}$ for this scenario.

For block diagonal covariance scenario, we used Gaussian $\mathbf{X}$
where $E(X_{im}) = 0$ and $Var(X_{im}) = 2$ for all $m$ with a block
diagonal covariance matrix. The $\mathbf{X}_m$’s were split into $G$
groups where for $m \neq k$ $Cov(\mathbf{X}_m,\mathbf{X}_k)=1$ if
predictor $m$ and $k$ were in the same group, otherwise
$Cov(\mathbf{X}_m,\mathbf{X}_k)=0$. Predictors with $\gamma_m=1$ were
randomly chosen from $G_1$ of the $G$ groups.

The methods use are:

- PROBE (all-at-once and one-at-a-time)

- Penalization methods:

  - LASSO (Tibshirani 1996),

  - Adaptive LASSO (Zou and Hastie 2005),

  - SCAD (Fan and Li 2001),

  - MCP (Zhang et al. 2010),

- SSLASSO (Ročková and George 2018)

- varbvs (Carbonetto and Stephens 2012)

- sparsevb (Ray and Szabó 2022)

- ebreg (Martin, Mess, and Walker 2017).

To fit the penalization methods, the penalty parameter that minimized
the $10$-fold cross-validated (CV) MSE was found using `glmnet`
(Friedman, Hastie, and Tibshirani 2010) for LASSO and ALASSO, and
`ncvreg` (Breheny and Huang 2011) for SCAD and MCP. The following
settings were used to fit the Bayesian variables selection methods (the
notation corresponds to the package documentation). To fit EMVS we used
a sequence from $-2$ to $2$ with increments of $0.25$ for $\log(\nu_0)$,
$\nu_1=1000$, a fixed prior on $\theta$ with $\theta=M_1/M$ (the true
value), and a conjugate prior on $\mathbf{\beta}$ and $\sigma^2$
(results correspond to $\nu_0=e^{-2}$). For SSLASSO we used
$\lambda_1=0.1$, $\lambda_0$ a series between $\lambda_1$ and $M_1$ with
400 elements, $a=M_1$ and $b=M-M_1$ with unknown variance. For VARBVS
the residual variance and prior variance of the regression coefficients
were approximated using maximum-likelihood. For SPARSEVB, Laplace and
Gaussian prior distributions were fit (results for Laplace were
marginally better and were used), the true value of $\sigma$ was used
for the residual standard deviation, and an intercept was included. For
EBREG we used $\alpha=0.99$, $\gamma=0.005$, a prior for $\sigma^2$ and
$\log\{I(|\mathbf{\gamma}| \leq n)/n\}$ as the log-prior on the model
size.

The simulation functions will automatically create and output a `csv`
file with the full results. Additionally, the function will summarized
results by iteration (if `verbose=TRUE`, the default) and print them in
the console.

## Squared Exponential Results

The simulations were ran for all 54 combinations of parameters given in
the arg_list.csv file for binary and continuous $X$. See the function
“simulation_func” in the “Sim_funcs.R” file for the code that runs the
simulation. Note that the `ebreg` function takes markedly longer to run
than the other functions. To run the simulation without ebreg, use
`ebreg_I = FALSE` when running simulation_func. In the paper, all
methods were ran for 1000 iterations except `ebreg`, which was ran for
100 iterations. Also, `ebreg` was not ran for M=400 due to the program
frequently returning an error in those settings.

Note: the simulations in the paper used the RandomFields (version
3.3.14) and RandomFieldsUtils (1.2.5) togenerate the MVN data with
squared exponential covariance structure. These packages have since been
removed from CRAN. For Windows machines, the corresponding versions of
RandomFields and RandomFieldsUtils can be installed from source files
from the archive at:

- <https://cran.r-project.org/web/packages/RandomFields/index.html>

- <https://cran.r-project.org/web/packages/RandomFieldsUtils/index.html>

For macOS, RandomFields (and all dependencies) can be installed via
macport:

- <https://ports.macports.org/port/R-RandomFields/>

`geoR` is an alternative package that has the same capabilities and can
be ran with the same parameters, however, it is **much** slower. For
`M=2500` generation of the data takes about 3 minutes and for `M=10000`
`M=10000` generation of the data takes multiple hours (with
`RandomFields` it takes around a minute). The simulation code will run
when either `RandomFields` or `geoR` are available.

## Block Diagonal Results

The simulations were ran for all 30 combinations of parameters given in
the arg_list_block.csv. These simulations were ran for all penalization
methods, SPARSEVB, and SSLASSO only.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-BreHua11" class="csl-entry">

Breheny, Patrick, and Jian Huang. 2011. “Coordinate Descent Algorithms
for Nonconvex Penalized Regression, with Applications to Biological
Feature Selection.” *Annals of Applied Statistics* 5 (1): 232–53.

</div>

<div id="ref-CarSte12" class="csl-entry">

Carbonetto, Peter, and Matthew Stephens. 2012. “Scalable Variational
Inference for Bayesian Variable Selection in Regression, and Its
Accuracy in Genetic Association Studies.” *Bayesian Analysis* 7 (1):
73–108. <https://doi.org/10.1214/12-BA703>.

</div>

<div id="ref-FanLi01" class="csl-entry">

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave
Penalized Likelihood and Its Oracle Properties.” *Journal of the
American Statistical Association* 96 (456): 1348–60.

</div>

<div id="ref-Frietal10" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software* 33 (1): 1–22.
<http://www.jstatsoft.org/v33/i01/>.

</div>

<div id="ref-Maretal17" class="csl-entry">

Martin, Ryan, Raymond Mess, and Stephen G Walker. 2017. “Empirical Bayes
Posterior Concentration in Sparse High-Dimensional Linear Models.”
*Bernoulli* 23 (3): 1822–47.

</div>

<div id="ref-RaySza22" class="csl-entry">

Ray, Kolyan, and Botond Szabó. 2022. “Variational Bayes for
High-Dimensional Linear Regression with Sparse Priors.” *Journal of the
American Statistical Association* 117 (539): 1270–81.

</div>

<div id="ref-RocGeo18" class="csl-entry">

Ročková, Veronika, and Edward I George. 2018. “The Spike-and-Slab
Lasso.” *Journal of the American Statistical Association* 113 (521):
431–44.

</div>

<div id="ref-Tib96" class="csl-entry">

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the
Lasso.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 58 (1): 267–88.

</div>

<div id="ref-Zha10" class="csl-entry">

Zhang, Cun-Hui et al. 2010. “Nearly Unbiased Variable Selection Under
Minimax Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942.

</div>

<div id="ref-ZouHas05" class="csl-entry">

Zou, Hui, and Trevor Hastie. 2005. “Regularization and Variable
Selection via the Elastic Net.” *Journal of the Royal Statistical
Society: Series B (Statistical Methodology)* 67 (2): 301–20.

</div>

</div>
