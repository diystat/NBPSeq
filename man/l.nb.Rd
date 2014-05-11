\name{l.nb}
\alias{l.nb}
\title{(private) The Log Likelihood of a NB Model}
\usage{
l.nb(kappa, mu, y)
}
\arguments{
  \item{kappa}{shape/size parameter}

  \item{mu}{mean parameter}

  \item{y}{a n-vector of NB counts}
}
\value{
the log likelihood of the NB model parameterized by
\code{(kappa, mu)}
}
\description{
The log likelihood of the NB model under the mean shape
parameterization
}
\details{
This function call dnbinom to compute the log likelihood
from each data point and sum the results over all data
points. kappa, mu and y should have compatible dimensions.
}

