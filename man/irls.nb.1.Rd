\name{irls.nb.1}
\alias{irls.nb.1}
\title{Estimate the regression coefficients in an NB GLM model}
\usage{
irls.nb.1(y, s, x, phi, beta0 = rep(NA, p), mustart = NULL, maxit = 50,
  tol.mu = 0.001/length(y), print.level = 0)
}
\arguments{
  \item{y}{an n vector of counts}

  \item{s}{a scalar or an n vector of effective library
  sizes}

  \item{x}{an n by p design matrix}

  \item{phi}{a scalar or an n-vector of dispersion
  parameters}

  \item{beta0}{a vector specifying known and unknown
  components of the regression coefficients: non-NA
  components are hypothesized values of beta, NA components
  are free components}

  \item{mustart}{starting values for the vector of means}

  \item{maxit}{maximum number of iterations}

  \item{tol.mu}{a number, convergence criteria}

  \item{tol}{a number, will be passed to Cdqrls}

  \item{print.level}{a number, print level}
}
\value{
a list of the following components: \item{beta}{a p-vector
of estimated regression coefficients} \item{mu}{an n-vector
of estimated mean values} \item{conv}{logical. Was the IRLS
algorithm judged to have converged?} \item{zero}{logical.
Was any of the fitted mean close to 0?}
}
\description{
Estimate the regression coefficients in an NB GLM model
with known dispersion parameters
}
\details{
This function estimates the regression coefficients using
iterative reweighted least squares (IRLS) algorithm, which
is equivalent to Fisher scoring. The implementation is
based on \code{glm.fit}.

Users can choose to fix some regression coefficients by
specifying \code{beta0}. (This is useful when fitting a
model under a null hypothesis.)
}

