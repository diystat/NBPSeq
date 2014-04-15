\name{estimate.norm.factors}
\alias{estimate.norm.factors}
\title{Estiamte Normalization Factors}
\usage{
estimate.norm.factors(counts, lib.sizes = colSums(counts),
  method = "AH2010")
}
\arguments{
  \item{counts}{a matrix of RNA-Seq read counts with rows
  corresponding to gene features and columns corresponding
  to independent biological samples.}

  \item{lib.sizes}{a vector of observed library sizes,
  usually and by default estimated by column totals.}

  \item{method}{a character string specifying the method
  for normalization, currenlty, can be NULL or "AH2010". If
  method=NULL, the normalization factors will have values
  of 1 (i.e., no normalization is applied); if
  method="AH2010" (default), the normalization method
  proposed by Anders and Huber (2010) will be used.}
}
\value{
a vector of normalization factors.
}
\description{
\code{estimate.norm.factors} estiamtes normalization
factors to account for apparent reduction or increase in
relative frequencies of non-differentially expressing genes
as a result of compensating the increased or decreased
relative frequencies of truly differentially expressing
genes.
}
\details{
We take gene expression to be indicated by relative
frequency of RNA-Seq reads mapped to a gene, relative to
library sizes (column sums of the count matrix). Since the
relative frequencies sum to 1 in each library (one column
of the count matrix), the increased relative frequencies of
truly over expressed genes in each column must be
accompanied by decreased relative frequencies of other
genes, even when those others do not truly differentially
express. If not accounted for, this may give a false
impression of biological relevance (see, e.g., Robinson and
Oshlack (2010), for some examples.)  A simple fix is to
compute the relative frequencies relative to effective
library sizes---library sizes multiplied by normalization
factors.
}
\examples{
## Load Arabidopsis data
data(arab)

## Estimate normalization factors using the method of Anders and Huber (2010)
norm.factors = estimate.norm.factors(arab);
print(norm.factors);
}
\references{
{Anders, S. and W. Huber (2010): "Differential expression
analysis for sequence count data," Genome Biol., 11, R106.

Robinson, M. D. and A. Oshlack (2010): "A scaling
normalization method for differential expression analysis
of RNA-seq data," Genome Biol., 11, R25.}
}

