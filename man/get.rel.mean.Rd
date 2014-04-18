\name{get.rel.mean}
\alias{get.rel.mean}
\title{(private) Extract row relative means of the pseudo counts for the specified group from an \code{nbp} object.}
\usage{
get.rel.mean(obj, grp.id)
}
\arguments{
  \item{obj}{a list with class \code{nbp}, output
  \code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
  \code{\link{exact.nb.test}} or \code{\link{nbp.test}}}

  \item{grp.id}{a number or a charater (same type as
  \code{obj$grp.ids}), group id}
}
\description{
(private) Extract row relative means of the pseudo counts
for the specified group from an \code{nbp} object.
}

