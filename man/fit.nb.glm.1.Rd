\name{fit.nb.glm.1}
\alias{fit.nb.glm.1}
\title{Fit a single negative binomial (NB) log-linear regression model
with known dispersion paramreters}
\usage{
fit.nb.glm.1(y, s, x, phi, beta0 = rep(NA, dim(x)[2]), ...)
}
\arguments{
  \item{y}{an n-vector of NB counts.}

  \item{s}{an n-vector of library sizes (multiplicative
  offset).}

  \item{x}{an n by p design matrix.}

  \item{phi}{a scalar or an n-vector, the NB dipsersion
  parameter.}

  \item{beta0}{a p-vector specifying the known and unknown
  components of beta, the regression coefficients. NA
  values indicate unknown components and non-NA values
  specify the values of the known components. The default
  is that all components of beta are unknown.}

  \item{...}{furhter arguements to be passed to
  \code{\link{irls.nb.1}}.}
}
\value{
a list \item{mu}{an n-vector, estimated means (MLE).}
\item{beta}{an p-vector, estimated regression coefficients
(MLE).} \item{iter}{number of iterations performed in the
IRLS algorithm.} \item{zero}{logical, whether any of the
estimated \code{mu} is close to zero.} \item{l}{log
likelihood of the fitted model.} \item{D}{a p-vector, the
score vector} \item{i}{a p-by-p matrix, fisher information
matrix} \item{j}{a p-by-p matrix, observed information
matrix}
}
\description{
Fit a NB log-linear regression model: find the MLE of the
regression coefficients and compute likelihood of the
fitted model, the score vector, and the Fisher and observed
information.
}
\details{
Under the NB regression model, the components of y follow a
NB distribution with means mu = s exp(x' beta) and
dispersion parameters phi.

The function will call \code{\link{irls.nb.1}} to find MLE
of the regression coefficients (which uses the iteratively
reweighted least squres (ILRS) algorithm).
}
\note{
The information matries, i and j, will be computed for all
all components of beta---including known components.
}

