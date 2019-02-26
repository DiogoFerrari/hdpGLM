hdpGLM
======

[![Travis build status](https://travis-ci.org/DiogoFerrari/hdpGLM.svg?branch=master)](https://travis-ci.org/DiogoFerrari/hdpGLM)

Overview
========

The package implements the hierarchical Dirichlet process Generalized
Linear Models proposed in the paper Modeling Context-Dependent Latent
Effect Heterogeneity. The model can be used to estimate latent
heterogeneity in the marginal effect of GLM linear coefficients, cluster
data points based on that latent heterogeneity, and investigate if
Simpson’s Paradox occurs due to latent or omitted features. It also can
be used with hierarchical data to estimate the effect of upper-level
covariates on the latent heterogeneity in the effect of lower-level
features.

For details of the model and the MCMC algorithm, see my paper [Modeling
Context-Dependent Latent Effect
Heterogeneity](https://dioferrari.files.wordpress.com/2018/09/hdpglm_v31.pdf),
current R&R at [Political
Analysis](https://www.cambridge.org/core/journals/political-analysis).

Instalation
===========

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
devtools::install_github("DiogoFerrari/hdpGLM")
# If you don't want to update the dependencies, use: (you may need to install some dependencies manually)
devtools::install_github("DiogoFerrari/hdpGLM", dependencies=F)
```

NOTE: it may be necessary to create a token to install the package from
the git repository in case it is private (see note at the bottom of help
page in R by running `help(install_github)`).

Usage
=====

Here is a simple example (for more information, see `help(hdpGLM)`).

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}

set.seed(10)
K    = 3 # number of latent clusters
nCov = 3 # number of observed covariates
simdata = hdpGLM_simulateData(400, nCov=nCov, K=K, family='gaussian')
data    = simdata$data
mcmc    = list(burn.in=1, n.iter=400)
samples = hdpGLM(y~., data=data, mcmc=mcmc, n.display=200)

summary(samples)
plot(samples)
plot(samples, terms="X1")
plot(samples, separate=T)
plot(samples, true.beta=summary(simdata)$beta)
plot(samples, true.beta=summary(simdata)$beta, separate=T)

```
