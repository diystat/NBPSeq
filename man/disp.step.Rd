\name{disp.step}
\alias{disp.step}
\title{(private) Specify a piecewise constant dispersion model}
\usage{
disp.step(counts, eff.lib.sizes, x, df = 1, knots = NULL, phi.pre = 0.1,
  mu.lower = 1, mu.upper = Inf)
}
\arguments{
  \item{counts}{a mxn matrix of NB counts}

  \item{eff.lib.sizes}{a n-vector of estimated effective
  library sizes}

  \item{x}{a nxp matrix, design matrix (specifying the
  treatment structure)}

  \item{df}{the number of steps}

  \item{knots}{a numerical vector of length df-1, giving
  the knots or jump locations.}

  \item{phi.pre}{a number, a preliminary constant
  dispersion value}

  \item{mu.lower}{a number, rows with mu.pre < mu.lower
  will not be used for estimating the dispersion model}

  \item{mu.upper}{a number, rows with mu.pre > mu.upper
  will not be used for estimating the dispersion model}
}
\value{
a list \item{fun}{a function that takes par as input and
outputs a matrix (same dimension as \code{counts}) of
dispersion values} \item{par.init}{a vector of length
\code{df}, initial values of par} \item{subset}{a logical
vector of length m, specifying a subset of genes to be used
when estimating the model parameters. Note that the
estimated model will be applied to all rows whenever
possible, but only rows specified in \code{subset} will be
used to estimate the dispersion model parameters.}
\item{pi.pre}{a m by n matrix, preliminary estimates of the
relative frequencies.} \item{knots}{a vector, the break
points of the step function}
}
\description{
Specify a piecewise constant (step) dispersion model. The
specified model are to be estimated using the function
\code{optim.disp.pl} or \code{optim.disp.apl}.
}
\details{
Under this model, the dispersion is modeled as a step
(piecewise constant) function.
}

