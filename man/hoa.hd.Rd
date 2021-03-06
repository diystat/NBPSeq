% Generated by roxygen2 (4.0.0): do not edit by hand
\name{hoa.hd}
\alias{hoa.hd}
\title{(private) HOA test for regression coefficients in an NBP GLM model}
\usage{
hoa.hd(y, s, x, phi, beta0, tol.mu = 0.001/length(y), print.level = 1)
}
\arguments{
\item{y}{an n vector of counts}

\item{s}{an n vector of effective library sizes}

\item{x}{an n by p design matrix}

\item{phi}{an n vector of dispersion parameters}

\item{beta0}{a p vector specifying null hypothesis: non NA
components are hypothesized values of beta, NA components are free
components}

\item{tol.mu}{convergence criteria}

\item{print.level}{a number, print level}
}
\value{
test statistics and  p-values of HOA, LR, and Wald tests
}
\description{
(private) HOA test for regression coefficients in an NBP GLM model
}

