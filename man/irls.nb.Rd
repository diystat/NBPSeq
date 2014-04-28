\name{irls.nb}
\alias{irls.nb}
\title{(private) Estiamte the regression coefficients in an NB GLM model}
\usage{
irls.nb(y, s, x, phi, beta0 = rep(NA, ncol(x)), mustart = NULL, ...,
  print.level = 0)
}
\arguments{
  \item{y}{an m*n matrix of counts}

  \item{s}{an n vector of effective library sizes}

  \item{x}{an n*p design matrix}

  \item{phi}{a scalar or an m*n matrix of NB2 dispersion
  coefficients}

  \item{beta0}{a K vector, non NA components are
  hypothesized values of beta, NA components are free
  components}

  \item{mustart}{an m*n matrix of starting values of the
  means}

  \item{...}{other parameters}

  \item{print.level}{a number, print level}
}
\value{
beta a K vector, the MLE of the regression coefficients.
}
\description{
Estimate the regression coefficients in an NBP GLM model
for each gene
}

