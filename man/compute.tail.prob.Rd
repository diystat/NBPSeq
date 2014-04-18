\name{compute.tail.prob}
\alias{compute.tail.prob}
\title{(private) Compute the tail probability of a conditional distribution
involving a pair of Negative Binomial (NB) random variables given
their sum}
\usage{
compute.tail.prob(s1, s2, mu1, mu2, kappa1, kappa2)
}
\arguments{
  \item{s1}{a number, the observed value of a NB random
  variable}

  \item{s2}{a number, the observed value of a NB random
  variable}

  \item{mu1}{a number, the mean parameter of the NB
  variable s1}

  \item{mu2}{a number, the mean parameter of the NB
  variable s2}

  \item{kappa1}{a number, the size parameter of the NB
  variable s1}

  \item{kappa2}{a number, the size parameter of the NB
  variable s2}
}
\value{
a number giving the probability of observing a (S1, S2)
that is as or more extreme than (s1, s2) given that
S1+S2=s1+s2.
}
\description{
Compute the probability of observing values of (S1, S2)
that are more extreme than (s1, s2) given that S1+S2=s1+s2
for a pair of Negative Binomial (NB) random variables (S1,
S2) with mean and size parameters (mu1, kappa1) and (mu2,
kappa2) respectively.
}
\details{
This function computes the probabily of (S1, S2) for all
values of S1 and S2 such that S1+S2=s1+s2, then sums over
the probabilites that are less than or equal to that of the
observed values (s1, s2). In context of DE test using
RNA-Seq data after thinning, S1 and S2 are often sums of
iid NB random variables (and are thus NB random variables
too).

The current implementation can be slow if s1 + s2 is large.

Potential improvements: For computing the one-sided tail
probability of Pr(S1<s1 |S1+S2=s1+s2), there might be a
faster way. The conditional distribution can be also
approximated by saddlepoint methods. If S1 and S2 are sum
of two subsets of iid random variables, the saddle point
approximation would be very accurate.
}

