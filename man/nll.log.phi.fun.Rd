\name{nll.log.phi.fun}
\alias{nll.log.phi.fun}
\title{Negative profile conditional likelihood of the dispersion model}
\usage{
nll.log.phi.fun(par, log.phi.fun, y, ls, n.grps, grps, grp.sizes, mu.lower,
  mu.upper, print.level)
}
\arguments{
  \item{par}{parameter of the disperesion model}

  \item{log.phi.fun}{the disperison model}

  \item{y}{counts}

  \item{ls}{library sizes}

  \item{n.grps}{number of groups}

  \item{grps}{a boolean matrix of group membership}

  \item{grp.sizes}{group sizes}

  \item{mu.lower}{lower bound for mu}

  \item{mu.upper}{upper bound for mu}

  \item{print.level}{print level}
}
\value{
negative log likelihood of the dispersion function
}
\description{
Negative profile conditional likelihood of the dispersion
model
}

