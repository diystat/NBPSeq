\name{fit.nb.glm.1u}
\alias{fit.nb.glm.1u}
\title{Fit a single negative binomial (NB) log-linear regression model with a common unknown dispersion paramreter}
\usage{
fit.nb.glm.1u(y, s, x, phi = NA, beta0 = rep(NA, dim(x)[2]),
  kappa = 1/phi, info.kappa = TRUE, ...)
}
\arguments{
  \item{y}{a n-vector of NB counts.}

  \item{s}{a n-vector of library sizes.}

  \item{x}{a n by p design matrix.}

  \item{phi}{a scalar, the NB dipsersion parameter.}

  \item{beta0}{a p-vector specifying the known and unknown
  components of beta, the regression coefficients. NA
  values indicate unknown components and non-NA values
  specify the values of the known components. The default
  is that all components of beta are unknown.}

  \item{kappa}{a scalar, the size/shape parameter.
  \code{kappa} will be set to \code{1/phi} if \code{phi} is
  not \code{NA} and will be estiamted if both \code{phi}
  and \code{kappa} are NA.}

  \item{info.kappa}{}

  \item{...}{additional parameters to
  \code{\link{irls.nb.1}}}
}
\value{
a list \item{mu}{an n-vector, estimated means (MLE).}
\item{beta}{an p-vector, estimated regression coefficients
(MLE).} \item{iter}{number of iterations performed in the
IRLS algorithm.} \item{zero}{logical, whether any of the
estimated \code{mu} is close to zero.} \item{kappa}{a
scalar, the size parameter} \item{phi}{a scalr, 1/kappa,
the dispsersion parameter} \item{l}{log likelihood of the
fitted model.} \item{D}{a p-vector, the score vector}
\item{j}{a p-by-p matrix, observed information matrix}
}
\description{
Fit a single negative binomial (NB) log-linear regression
model with a common unknown dispersion paramreter.
}
\details{
Find the MLE of the dipsersion parameter and the regression
coefficients in a NB regression model.

Under the NB regression model, the components of y follow a
NB distribution with means mu = s exp(x' beta) and a common
dispersion parameter phi.
}
\note{
When the disperison is known, the user should specify only
one of \code{phi} or \code{kappa}. Whenever \code{phi} is
specified (non-NA), \code{kappa} will be set to
1/\code{phi}.

The observed information matrix, j, will be computed for
all parameters---kappa and all components of beta
(including known components). It will be computed at the
estimated values of (phi, beta) or (kappa, beta), which can
be unconstrained or constrained MLEs depending on how the
arguments \code{phi} (or \code{kappa}) and \code{beta} are
specified.

TODO: allow computing the information matrix using phi or
log(kappa) as parameter
}

