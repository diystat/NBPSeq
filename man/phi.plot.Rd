\name{phi.plot}
\alias{phi.plot}
\title{Plot estimated genewise NB2 dispersion parameter versus
estimated mean}
\usage{
phi.plot(counts, alpha = 2, xlab = "mean", ylab = "phi.hat",
  main = "phi.hat vs mean", log = "xy", ...)
}
\arguments{
  \item{counts}{a matrix of NB counts}

  \item{alpha}{alpha}

  \item{xlab}{x label}

  \item{ylab}{y label}

  \item{main}{main}

  \item{log}{log}

  \item{...}{other}
}
\description{
Plot estimated NB2 dispersion parameter versus estimated
mean
}
\details{
\code{phi.plot} estimate the NB2 dispersion parameter for
each gene separately by \eqn{\phi = (v - \mu) / \mu^alpha},
where \eqn{\mu} and \eqn{v} are sample mean and sample
variance. By default, \eqn{alpha=2}.
}
\note{
Currently, we discards genes giving 0 mean or negative
dispersion estimate (which can happen if sample variance is
smaller than the sample mean).
}

