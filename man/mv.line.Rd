\name{mv.line}
\alias{mv.line}
\title{(private) Overlay an estimated mean-variance line}
\usage{
mv.line(mu, v, ...)
}
\arguments{
  \item{mu}{a vector of mean values}

  \item{v}{a vector of variance values}

  \item{...}{other}
}
\description{
Overlay an estimated mean-variance line on an existing
mean-variance plot
}
\details{
Users should call mv.plot before calling this function.

If the length of theinput vectors (\code{mu}, \code{v}) is
greater than 1000, then we will only use a subset of the
input vectors.
}

