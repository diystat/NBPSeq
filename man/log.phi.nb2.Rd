\name{log.phi.nb2}
\alias{log.phi.nb2}
\title{(private) A NBP dispersion model}
\usage{
\method{log}{phi.nb2}(par, pi)
}
\arguments{
  \item{par}{a number, log dispersion}
}
\value{
a vector of length m, log dipserison.
}
\description{
Specify a NB2 dispersion model. The parameter of the
specified model are to be estimated from the data using the
function \code{optim.pcl}.
}
\details{
Under this NB2 model, the log dispersion is a constant.
}

