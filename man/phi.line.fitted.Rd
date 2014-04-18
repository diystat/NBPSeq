\name{phi.line.fitted}
\alias{phi.line.fitted}
\title{(private) Overlay an estimated mean-dispersion line on an existing plot}
\usage{
phi.line.fitted(obj, alpha = 2, ...)
}
\arguments{
  \item{obj}{a list with two components: \code{mu}, a
  vector of mean values; \code{v}, a vector of variance
  values.}

  \item{alpha}{alpha}

  \item{...}{other}
}
\description{
Overlay an estimated mean-dispersion line on an existing
plot
}
\details{
This function is a wrapper of \code{\link{phi.line}}. It
takes a list (rather than two separate vectors) as input.
}
\note{
Users should call phi.plot before calling this function.
}
\seealso{
\code{\link{phi.line}}
}

