\name{optim.disp.pl}
\alias{optim.disp.pl}
\title{(private) Estimate the parameters in a dispersion model}
\usage{
optim.disp.pl(disp, counts, eff.lib.sizes, x, method = "L-BFGS-B", ...)
}
\arguments{
  \item{disp}{a list, output from \code{\link{disp.nbp}},
  \code{\link{disp.nbq}}, \code{\link{disp.nbs}}, and so
  on.}

  \item{counts}{a matrix, the nb counts}

  \item{eff.lib.sizes}{effective library sizes}

  \item{x}{a desing matrix}

  \item{method}{the optimization method to be used by
  \code{optim}}

  \item{...}{additional parareters, will be passed to
  optim().}
}
\value{
a list with components:
}
\description{
Estimate the parameters in a dispersion model.
}
\details{
The function will call the R funciton \code{optim} to
mimimize the negative log likelihood of the dipserison
model.
}

