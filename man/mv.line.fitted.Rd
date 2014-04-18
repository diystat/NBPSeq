\name{mv.line.fitted}
\alias{mv.line.fitted}
\title{(private) Overlay an estimated mean-variance line}
\usage{
mv.line.fitted(obj, ...)
}
\arguments{
  \item{obj}{a list with components mu, a vector of mean
  values, and v, a vector of variance values.}

  \item{...}{other parameters}
}
\description{
Overlay an estimated mean-variance line on existing plot
}
\details{
This functions is a wrapper of \code{\link{mv.line}}. It
takes a list (rather than two vectors) as input.
}
\note{
Users should call mv.plot before calling this function.
}
\seealso{
\code{\link{mv.line}}
}

