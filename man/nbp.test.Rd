\name{nbp.test}
\alias{nbp.test}
\title{NBP Test for Differential Gene Expression from RNA-Seq
Counts}
\usage{
nbp.test(counts, grp.ids, grp1, grp2, norm.factors = rep(1, dim(counts)[2]),
  model.disp = "NBQ", lib.sizes = colSums(counts), print.level = 1, ...)
}
\arguments{
  \item{counts}{an \eqn{n} by \eqn{r} matrix of RNA-Seq
  read counts with rows corresponding to genes (exons, gene
  isoforms, etc) and columns corresponding to libraries
  (independent biological samples).}

  \item{grp.ids}{an \eqn{r} vector of treatment group
  identifiers (e.g. integers).}

  \item{grp1}{group 1 id}

  \item{grp2}{group 2 id}

  \item{norm.factors}{an \eqn{r} vector of normalization
  factors.}

  \item{model.disp}{a string, one of "NB2", "NBP" or "NBQ"
  (default).}

  \item{lib.sizes}{(unnormalized) library sizes}

  \item{print.level}{a number, controls the amount of
  messages printed: 0 for suppressing all messages, 1
  (default) for basic progress messages, and 2 to 5 for
  increasingly more detailed messages.}

  \item{...}{optional parameters to be passed to
  \code{\link{estimate.disp}}, the function that estimates
  the dispersion parameters.}
}
\value{
a list with the following components: \item{counts}{an
\eqn{n} by \eqn{r} matrix of counts, same as input.}
\item{lib.sizes}{an \eqn{r} vector, column sums of the
count matrix.} \item{grp.ids}{an \eqn{r} vector,
identifiers of treatment groups, same as input.}
\item{grp1, grp2}{identifiers of the two groups to be
compared, same as input.} \item{eff.lib.sizes}{an \eqn{r}
vector, effective library sizes, lib.sizes multiplied by
the normalization factors.} \item{pseudo.counts}{count
matrix after thinning, same dimension as counts}
\item{pseduo.lib.sizes}{an \eqn{r} vector, effective
library sizes of pseudo counts, i.e., column sums of the
pseudo count matrix multiplied by the normalization.}
\item{phi, alpha}{two numbers, parameters of the dispersion
model.} \item{pie}{a matrix, same dimension as
\code{counts}, estimated mean relative frequencies of
RNA-Seq reads mapped to each gene.} \item{pooled.pie}{a
matrix, same dimenions as \code{counts}, estimated pooled
mean of relative frequencies in the two groups being
compared.} \item{expression.levels}{a \eqn{n} by 3 matrix,
estimated gene expression levels as indicated by mean
relative frequencies of RNA-Seq reads. It has three columns
\code{grp1}, \code{grp2}, \code{pooled} corresponding to
the two treatment groups and the pooled mean.}
\item{log.fc}{an \eqn{n}-vector, base 2 log fold change in
mean relative frequency between two groups.}
\item{p.values}{an \eqn{n}-vector, p-values of the exact NB
test applied to each gene (row).} \item{q.values}{an
\eqn{n}-vector, q-values (estimated FDR) corresponding to
the p-values.}
}
\description{
\code{nbp.test} fits an NBP model to the RNA-Seq counts and
performs Robinson and Smyth's exact NB test on each gene to
assess differential gene expression between two groups.
}
\details{
\code{nbp.test} calls \code{\link{prepare.nbp}} to create
the NBP data structure, perform optional normalization and
adjust library sizes, calls \code{\link{estimate.disp}} to
estimate the NBP dispersion parameters and
\code{\link{exact.nb.test}} to perform the exact NB test
for differential gene expression on each gene. The results
are summarized using p-values and q-values (FDR).

\subsection{Overview}{ For assessing evidence for
differential gene expression from RNA-Seq read counts, it
is critical to adequately model the count variability
between independent biological replicates.  Negative
binomial (NB) distribution offers a more realistic model
for RNA-Seq count variability than Poisson distribution and
still permits an exact (non-asymptotic) test for comparing
two groups.

For each individual gene, an NB distribution uses a
dispersion parameter \eqn{\phi_i} to model the
extra-Poisson variation between biological replicates.
Across all genes, parameter \eqn{\phi_i} tends to vary with
the mean \eqn{\mu_i}. We capture the dispersion-mean
dependence using a parametric model: NB2, NBP and NBQ. (See
\code{\link{estimate.disp}} for more details.)}

\subsection{Count Normalization}{ We take gene expression
to be indicated by relative frequency of RNA-Seq reads
mapped to a gene, relative to library sizes (column sums of
the count matrix). Since the relative frequencies sum to 1
in each library (one column of the count matrix), the
increased relative frequencies of truly over expressed
genes in each column must be accompanied by decreased
relative frequencies of other genes, even when those others
do not truly differentially express. Robinson and Oshlack
(2010) presented examples where this problem is noticeable.

A simple fix is to compute the relative frequencies
relative to effective library sizes---library sizes
multiplied by normalization factors.  By default,
\code{nbp.test} assumes the normalization factors are 1
(i.e. no normalization is needed). Users can specify
normalization factors through the argument
\code{norm.factors}.  Many authors (Robinson and Oshlack
(2010), Anders and Huber (2010)) propose to estimate the
normalization factors based on the assumption that most
genes are NOT differentially expressed. }

\subsection{Library Size Adjustment}{ The exact test
requires that the effective library sizes (column sums of
the count matrix multiplied by normalization factors) are
approximately equal. By default, \code{nbp.test} will thin
(downsample) the counts to make the effective library sizes
equal. Thinning may lose statistical efficiency, but is
unlikely to introduce bias.}
}
\note{
Due to thinning (random downsampling of counts), two
identical calls to \code{nbp.test} may yield slightly
different results. A random number seed can be used to make
the results reproducible. The regression analysis method
implemented in \code{\link{nb.glm.test}} does not require
thinning and can also be used to compare expression in two
groups.

Advanced users can call
\code{\link{estimate.norm.factors}},
\code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
\code{\link{exact.nb.test}} directly to have more control
over modeling and testing.
}
\examples{
## Load Arabidopsis data
data(arab);

## Specify treatment groups and ids of the two groups to be compared
grp.ids = c(1, 1, 1, 2, 2, 2);
grp1 = 1;
grp2 = 2;

## Estimate normalization factors
norm.factors = estimate.norm.factors(arab);

## Set a random number seed to make results reproducible
set.seed(999);

## Fit the NBP model and perform exact NB test for differential gene expression.
## For demonstration purpose, we will use the first 100 rows of the arab data.
res = nbp.test(arab[1:100,], grp.ids, grp1, grp2,
  lib.sizes = colSums(arab), norm.factors = norm.factors, print.level=3);

## The argument lib.sizes is needed since we only use a subset of
## rows. If all rows are used, the following will be adequate:
##
## res = nbp.test(arab, grp.ids, grp1, grp2, norm.factors = norm.factors);

## Show top ten most differentially expressed genes
subset = order(res$p.values)[1:10];
print(res, subset);

## Count the number of differentially expressed genes (e.g. qvalue < 0.05)
alpha = 0.05;
sig.res = res$q.values < alpha;
table(sig.res);

## Show boxplots, MA-plot, mean-variance plot and mean-dispersion plot
par(mfrow=c(3,2));
plot(res);
}
\references{
Di, Y, D. W. Schafer, J. S. Cumbie, and J. H. Chang (2011):
"The NBP Negative Binomial Model for Assessing Differential
Gene Expression from RNA-Seq", Statistical Applications in
Genetics and Molecular Biology, 10 (1).

Robinson, M. D. and G. K. Smyth (2007): "Moderated
statistical tests for assessing differences in tag
abundance," Bioinformatics, 23, 2881-2887.

Robinson, M. D. and G. K. Smyth (2008): "Small-sample
estimation of negative binomial dispersion, with
applications to SAGE data," Biostatistics, 9, 321-332.

Anders, S. and W. Huber (2010): "Differential expression
analysis for sequence count data," Genome Biol., 11, R106.

Robinson, M. D. and A. Oshlack (2010): "A scaling
normalization method for differential expression analysis
of RNA-seq data," Genome Biol., 11, R25.
}
\seealso{
\code{\link{prepare.nbp}}, \code{\link{estimate.disp}},
\code{\link{exact.nb.test}}.
}

