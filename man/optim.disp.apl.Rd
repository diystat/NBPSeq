\name{optim.disp.apl}
\alias{optim.disp.apl}
\title{(private) Estimate the parameters in a dispersion model}
\usage{
optim.disp.apl(disp, counts, eff.lib.sizes, x, method = "L-BFGS-B",
  print.level = 1, ...)
}
\arguments{
  \item{disp}{a list, output from \code{\link{disp.nbp}},
  \code{\link{disp.nbq}}.}

  \item{counts}{a matrix, the nb counts}

  \item{eff.lib.sizes}{effective library sizes}

  \item{x}{a desing matrix}

  \item{method}{the optimization method to be used by
  \code{optim}}

  \item{print.level}{print level}

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
mimimize the negative log adjusted profile likelihood of
the dipserison model.
}

