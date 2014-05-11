##' Fit a NB log-linear regression model: find the MLE of the
##' regression coefficients and compute likelihood of the fitted model,
##' the score vector, and the Fisher and observed information.
##'
##' Under the NB regression model, the components of y follow a NB
##' distribution with means mu = s exp(x' beta) and dispersion
##' parameters phi.
##'
##' The function will call \code{\link{irls.nb.1}} to find MLE of the
##' regression coefficients (which uses the iteratively reweighted
##' least squres (ILRS) algorithm).
##'
##' @note
##'
##' The information matries, i and j, will be computed for all all components
##' of beta---including known components.
##'
##' @title Fit a single negative binomial (NB) log-linear regression model
##' with known dispersion paramreters
##'
##' @param y an n-vector of NB counts.
##' @param s an n-vector of library sizes (multiplicative offset).
##' @param x an n by p design matrix.
##' @param phi a scalar or an n-vector, the NB dipsersion parameter.
##' @param beta0 a p-vector specifying the known and unknown
##' components of beta, the regression coefficients. NA values
##' indicate unknown components and non-NA values specify the values
##' of the known components. The default is that all components of
##' beta are unknown.
##' @param ... furhter arguements to be passed to \code{\link{irls.nb.1}}.
##' @return a list
##'   \item{mu}{an n-vector, estimated means (MLE).}
##'   \item{beta}{an p-vector, estimated regression coefficients (MLE).}
##'   \item{iter}{number of iterations performed in the IRLS algorithm.}
##'   \item{zero}{logical, whether any of the estimated \code{mu} is close to zero.}
##'   \item{l}{log likelihood of the fitted model.}
##'   \item{D}{a p-vector, the score vector}
##'   \item{i}{a p-by-p matrix, fisher information matrix} 
##'   \item{j}{a p-by-p matrix, observed information matrix}
fit.nb.glm.1 = function(y, s, x, phi, beta0=rep(NA, dim(x)[2]), 
  ...) {

  ## Find MLE of beta
  res = irls.nb.1(y, s, x, phi, beta0, ...);

  kappa = 1/phi;
  mu = res$mu;

  ##
  res$phi =phi;
  
  ## Compute the log likelihod
  res$l = ll.nb(kappa, mu, y);

  v = drop(mu + phi * mu^2);

  ## Compute the score vector;
  n = length(y);
  res$D1 = matrix((y - mu) * mu/v, 1, n) %*% x;
  
  ## Compute Fisher and observed info
  res$i = t(x) %*% diag(mu^2/v) %*% x;
  res$j = t(x) %*% diag(mu^2 * (y/mu^2 - (y+kappa)/(mu+kappa)^2) - (y-mu)*mu/v) %*% x;

  res
}

##' Fit a single negative binomial (NB) log-linear regression model with a common unknown dispersion paramreter.
##'
##' Find the MLE of the dipsersion parameter and the regression
##' coefficients in a NB regression model.
##'
##' Under the NB regression model, the components of y follow a NB
##' distribution with means mu = s exp(x' beta) and a common
##' dispersion parameter phi.
##'
##' @note
##'
##' When the disperison is known, the user should specify only one of
##' \code{phi} or \code{kappa}. Whenever \code{phi} is specified
##' (non-NA), \code{kappa} will be set to 1/\code{phi}.
##'
##' The observed information matrix, j, will be computed for all
##' parameters---kappa and all components of beta (including known
##' components). It will be computed at the estimated values of (phi,
##' beta) or (kappa, beta), which can be unconstrained or constrained
##' MLEs depending on how the arguments \code{phi} (or \code{kappa})
##' and \code{beta} are specified.
##'
##' TODO: allow computing the information matrix using phi or log(kappa) as
##' parameter
##'
##' @title Fit a single negative binomial (NB) log-linear regression model with a common unknown dispersion paramreter
##'
##' @param y a n-vector of NB counts.
##' @param s a n-vector of library sizes.
##' @param x a n by p design matrix.
##' @param phi a scalar, the NB dipsersion parameter.
##' @param beta0 a p-vector specifying the known and unknown
##' components of beta, the regression coefficients. NA values
##' indicate unknown components and non-NA values specify the values
##' of the known components. The default is that all components of
##' beta are unknown.
##' @param kappa a scalar, the size/shape parameter.  \code{kappa}
##' will be set to \code{1/phi} if \code{phi} is not \code{NA}  and
##' will be estiamted if both \code{phi} and \code{kappa} are NA.
##' @param info.kappa 
##' @param ... additional parameters to \code{\link{irls.nb.1}}
##' @return a list
##'   \item{mu}{an n-vector, estimated means (MLE).}
##'   \item{beta}{an p-vector, estimated regression coefficients (MLE).}
##'   \item{iter}{number of iterations performed in the IRLS algorithm.}
##'   \item{zero}{logical, whether any of the estimated \code{mu} is close to zero.}
##'   \item{kappa}{a scalar, the size parameter}
##'   \item{phi}{a scalr, 1/kappa, the dispsersion parameter}
##'   \item{l}{log likelihood of the fitted model.}
##'   \item{D}{a p-vector, the score vector}
##'   \item{j}{a p-by-p matrix, observed information matrix}
fit.nb.glm.1u = function(y, s, x, phi=NA, beta0=rep(NA, dim(x)[2]), kappa=1/phi,
  info.kappa=TRUE, ...) {

## TODO:  @param info.kappa logical, if \code{TRUE} the information matrix will be
## computed for (kappa, beta), otherwise the information matrix will
## be computed for (phi, beta).

  if (is.na(kappa)) {
    
    ## Find preliminary estiamtes of mu assuming phi=0.1. Will serve as
    ## initial values for the later irls algorithm.
    mustart = irls.nb.1(y, s, x, 0.1, beta0, ...)$mu;

    ## Log likelihood of log(kappa)
    ll =  function(lkappa) {
      kappa = exp(lkappa);
      res = irls.nb.1(y, s, x, 1/kappa, beta0, mustart);
      sum(dnbinom(y, size =kappa, mu=res$mu, log=TRUE)); 
      ## log.likelihood.nb(kappa, res$mu, y);
    }

    res = optimize(ll, c(log(1e-20), log(1e20)), maximum=TRUE);
    kappa.hat = exp(res$maximum);
  } else {
    kappa.hat = kappa
  }

  phi.hat = 1/kappa.hat;

  res = irls.nb.1(y, s, x, phi.hat, beta0, mustart, ...);

  beta.hat = res$beta;
  mu.hat = res$mu;

  ## Compute the observed information, j
  v.hat = drop(mu.hat + phi.hat * mu.hat^2);
  n = length(y);
  n.pars = length(beta0) + 1;

  j.hat = matrix(NA, n.pars, n.pars);
  j.hat[-1,-1] = t(x) %*% diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa.hat)/(mu.hat+kappa.hat)^2) - (y-mu.hat)*mu.hat/v.hat) %*% x;
  ## Be careful with the sign!
  j.hat[1,-1] = j.hat[-1,1] = - matrix((y - mu.hat) * mu.hat / (mu.hat+kappa.hat)^2, 1, n) %*% x;
  j.hat[1, 1] = -sum(trigamma(kappa.hat+y) - trigamma(kappa.hat) + 1/kappa.hat - 2/(kappa.hat+mu.hat) + (y+kappa.hat)/(mu.hat + kappa.hat)^2);

  ##  i.hat = t(x) %*% diag(mu.hat^2/v.hat) %*% x;
  ## res$l = ll.nb(kappa.hat, mu = mu.hat, y);

  ## See whether Poisson model fits better
  ## l.poisson  = sum(dpois(y, lambda = mu.hat, log=TRUE));
  ## eps = 1e-7;
  ## if (l.poisson + eps >= l & kappa.hat > 1e8) {
  ## l =  l.poisson;
  ##   kappa.hat = Inf;
  ##}
  ## l = sum(dnbinom(y, size = 1/phi, mu=res$mu, log=TRUE)); 

  res$phi = phi.hat
  res$kappa = kappa.hat;
  res$l = ll.nb(kappa.hat, mu.hat, y);
  res$j = j.hat

  res;
}

