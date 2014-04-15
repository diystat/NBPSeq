\name{prepare.nb.data}
\alias{prepare.nb.data}
\title{Prepare the NB Data Structure for RNA-Seq Read Counts}
\usage{
prepare.nb.data(counts, lib.sizes = colSums(counts), norm.factors = rep(1,
  dim(counts)[2]), tags = NULL)
}
\arguments{
  \item{counts}{an mxn matrix of RNA-Seq read counts with
  rows corresponding to gene features and columns
  corresponding to independent biological samples.}

  \item{lib.sizes}{an n-vector of observed library sizes.
  By default, library sizes are estimated to the column
  totals of the matrix \code{counts}.}

  \item{norm.factors}{an n-vector of normalization factors.
  By default, have values 1 (no normalization is applied).}

  \item{tags}{a matrix of tags associated with genes, one
  row for each gene (having the same number of rows as
  \code{counts}.}
}
\value{
A list containing the following components:
\item{counts}{the count matrix, same as input.}
\item{lib.sizes}{observed library sizes, same as input.}
\item{norm.factors}{normalization factors, same as input.}
\item{eff.lib.sizes}{effective library sizes
(\code{lib.sizes} x \code{norm.factors}).}
\item{rel.frequencies}{relative frequencies (counts divided
by the effective library sizes).} \item{tags}{a matrix of
gene tags, same as input.}
}
\description{
Create a data structure to hold the RNA-Seq read counts and
other relevant information.
}

