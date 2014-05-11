##' The log likelihood of the NB model under the mean shape parameterization
##'
##' This function call dnbinom to compute the log likelihood from each data point and sum the results over all data points.
##' kappa, mu and y should have compatible dimensions.
##' 
##' @title (private) The Log Likelihood of a NB Model
##'
##' @param kappa shape/size parameter
##' @param mu mean parameter
##' @param y a n-vector of NB counts
##' @return the log likelihood of the NB model parameterized by \code{(kappa, mu)}
l.nb= function(kappa, mu, y) {
  sum(dnbinom(y, kappa, mu=mu, log=TRUE));
}
