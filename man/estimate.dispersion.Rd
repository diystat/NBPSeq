% Generated by roxygen2 (4.0.0): do not edit by hand
\name{estimate.dispersion}
\alias{estimate.dispersion}
\title{Estimate Negative Binomial Dispersion}
\usage{
estimate.dispersion(nb.data, x, model = "NBQ", predictor = "pi",
  method = "MAPL", fast = TRUE, ...)
}
\arguments{
\item{nb.data}{output from \code{\link{prepare.nb.data}}.}

\item{x}{a design matrix specifying the mean structure of each row.}

\item{model}{the name of the dispersion model, one of "NB2", "NBP", "NBQ" (default), "NBS" or "step".}

\item{predictor}{}

\item{method}{a character string specifying the method for estimating the dispersion model, one of "ML" or "MAPL" (default).}

\item{fast}{use a faster (but might be less accurate method)}

\item{...}{additional parameters to optim.fun.<method>}
}
\value{
a list with following components:
\item{estimates}{dispersion estimates for each read count, a matrix of the same dimensions as
the \code{counts} matrix in \code{nb.data}.}
\item{likelihood}{the likelihood of the fitted model.}
\item{model}{details of the estimate dispersion model, NOT intended for use by end users. The name and contents of this component are subject to change in future versions.}
}
\description{
Estimate NB dispersion by modeling it as a parametric function of preliminarily estimated log mean relative frequencies.
}
\details{
We use a negative binomial (NB) distribution to model the read
frequency of gene \eqn{i} in sample \eqn{j}.  A negative binomial
(NB) distribution uses a dispersion parameter \eqn{\phi_{ij}} to
model the extra-Poisson variation between biological replicates.
Under the NB model, the mean-variance relationship of a single
read count satisfies \eqn{\sigma_{ij}^2 = \mu_{ij} + \phi_{ij}
\mu_{ij}^2}.  Due to the typically small sample sizes of RNA-Seq
experiments, estimating the NB dispersion \eqn{\phi_{ij}} for each
gene \eqn{i} separately is not reliable.  One can pool information
across genes and biological samples by modeling \eqn{\phi_{ij}} as
a function of the mean frequencies and library sizes.

Under the NB2 model, the dispersion is a constant across all genes and samples.

Under the NBP model, the log dispersion is modeled as a linear
function of the preliminary estimates of the log mean relative
frequencies (\code{pi.pre}):

log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset),

where \code{pi.offset} is 1e-4.

Under the NBQ model, the dispersion is modeled as a quadratic
function of the preliminary estimates of the log mean relative
frequencies (pi.pre):

log(phi) =  par[1] + par[2] * z + par[3] * z^2,

where z = log(pi.pre/pi.offset). By default, pi.offset is the median of pi.pre[subset,].

Under this NBS model, the dispersion is
modeled as a smooth function (a natural cubic spline function) of
the preliminary estimates of the log mean relative frequencies
(pi.pre).

Under the "step" model, the dispersion is modeled as a step
(piecewise constant) function.
}
\note{
Currently, it is unclear whether a dispersion-modeling
approach will outperform a more basic approach where regression
model is fitted to each gene separately without considering the
dispersion-mean dependence. Clarifying the power-robustness of the
dispersion-modeling approach is an ongoing research topic.
}
\examples{
## See the example for test.coefficient.
}

