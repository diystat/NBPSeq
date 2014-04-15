\name{print.nbp}
\alias{print.nbp}
\title{Print summary of an NBP Object}
\usage{
\method{print}{nbp}(x, subset = 1:10, ...)
}
\arguments{
  \item{x}{Output from \code{\link{prepare.nbp}},
  \code{\link{estimate.disp}}, or \code{\link{nbp.test}}.}

  \item{subset}{indices of rows of the count matrix to be
  printed.}

  \item{...}{other parameters (for future use).}
}
\description{
Print contents of an NBP object, output from
\code{\link{prepare.nbp}}, \code{\link{estimate.disp}}, or
\code{\link{nbp.test}}.
}
\examples{
## See the example for nbp.test
}
\seealso{
\code{\link{nbp.test}}.
}

