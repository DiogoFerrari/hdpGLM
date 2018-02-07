Overview
========

Hierarchical Dirichlet Process Generalized Linear models **hdpGLM** is a
package that provides an implementation of semi-parametric models,
particularly infinite mixture of generalized linear (GLM) or generalized
linear mixed (GLMM) models. Dirichlet Process priors are used both for
the mixture probabilities and, optionally, for the random effects. It is
a generalization of the classical GLM and GLMM. It is developed to
estimate clusters in the data generated by latent or unobserved
covariates when they exist. The model estimate clusters when there is
heterogeneity in the population in terms of how the vector of covariates
is associated with the outcome. When there is no heterogeneity, the
model results are similar to classical generalized linear model. The
package also provides a hierarchical version of the model in which the
number and characteristics of the clusters in different contexts can
depend on context level characteristics.

Instalation
===========

``` {.r .rundoc-block rundoc-language="R" rundoc-exports="code"}
devtools::install_github("hdpGLM")
```

Usage
=====

More information
================
