% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src_summaries.R
\name{plot.dpGLM}
\alias{plot.dpGLM}
\title{Default plot for class dpGLM}
\usage{
\method{plot}{dpGLM}(
  x,
  terms = NULL,
  separate = FALSE,
  hpd = TRUE,
  true.beta = NULL,
  title = NULL,
  subtitle = NULL,
  adjust = 1,
  ncols = NULL,
  only.occupied.clusters = TRUE,
  focus.hpd = FALSE,
  legend.position = "top",
  colour = "grey",
  alpha = 0.4,
  display.terms = TRUE,
  plot.mean = TRUE,
  legend.label.true.value = "True",
  ...
)
}
\arguments{
\item{x}{a dpGLM object with the samples from generated by \code{\link{hdpGLM}}}

\item{terms}{string vector with the name of covariates to plot. If \code{NULL} (default), all covariates are plotted.}

\item{separate}{boolean, if \code{TRUE} the linear coefficients \code{beta} will be displayed in their separate clusters.}

\item{hpd}{boolean, if \code{TRUE} and \code{separate=T}, the 95\% HPDI lines will be displayed.}

\item{true.beta}{either \code{NULL} (default) or a \code{data.frame} with the true values of the linear coefficients \code{beta} if they are known. The \code{data.frame} must contain a column named \code{k} indicating the cluster of \code{beta}, and a column named \code{Parameter} with the name of the linear coefficients (\code{beta1}, \code{beta2}, ..., \code{beta_dx}, where \code{dx} is the number of covariates at the individual level, and beta1 is the coefficient of the intercept term). It must contain a column named \code{True} with the true value of the \code{betas}.}

\item{title}{string, the title of the plot}

\item{subtitle}{string, the subtitle of the plot}

\item{adjust}{the bandwidth used is actually \code{adjust*bw}.
    This makes it easy to specify values like \sQuote{half the default}
    bandwidth.}

\item{ncols}{integer, the number of columns in the plot}

\item{only.occupied.clusters}{boolean, if \code{TRUE} it shows only the densities of the clusters that actually have data points assigned to it with high probability}

\item{focus.hpd}{boolean, if \code{TRUE} and separate is also \code{TRUE} it will display only the 95\% HPDI of the posterior density of the linear coefficients \code{beta}}

\item{legend.position}{one of four options: "bottom" (default), "top", "left", or "right". It indicates the position of the legend}

\item{colour}{= string with color to fill the density plot}

\item{alpha}{number between 0 and 1 indicating the degree of transparency of the density}

\item{display.terms}{boolean, if \code{TRUE} (default), the covariate name is displayed in the plot}

\item{plot.mean}{boolean, if \code{TRUE} the posterior mean of every cluster is displayed}

\item{legend.label.true.value}{a string with the value to display in the legend when the \code{true.beta} is used}

\item{...}{ignored}
}
\description{
This function generates desity plots with the posterior distribution generated by the function \code{\link{hdpGLM}}
}
\examples{
# Note: this example is just for illustration. MCMC iterations are very reduced
set.seed(10)
n = 20
data = tibble::tibble(x1 = rnorm(n, -3),
                                   x2 = rnorm(n,  3),
                                   z  = sample(1:3, n, replace=TRUE),
                                   y  =I(z==1) * (3 + 4*x1 - x2 + rnorm(n)) +
                                       I(z==2) * (3 + 2*x1 + x2 + rnorm(n)) +
                                       I(z==3) * (3 - 4*x1 - x2 + rnorm(n)) ,
                                   ) 

## estimation
mcmc    = list(burn.in=1, n.iter=50)
samples = hdpGLM(y ~ x1 + x2,  data=data, mcmc=mcmc, n.display=1)

plot(samples)


}
