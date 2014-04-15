\name{exact.nb.test}
\alias{exact.nb.test}
\title{Exact Negative Binomial Test for Differential Gene Expression}
\usage{
exact.nb.test(obj, grp1, grp2, print.level = 1)
}
\arguments{
  \item{obj}{output from \code{\link{estimate.disp}}.}

  \item{grp1}{Identifier of group 1. A number, character or
  string (should match at least one of the obj$grp.ids).}

  \item{grp2}{Identifier of group 2. A number, character or
  string (should match at least one of the obj$grp.ids).}

  \item{print.level}{a number.  Controls the amount of
  messages printed: 0 for suppressing all messages, 1 for
  basic progress messages, larger values for more detailed
  messages.}
}
\value{
the list \code{obj} from the input with the following added
components: \item{grp1}{same as input.} \item{grp2}{same as
input.} \item{pooled.pie}{estimated pooled mean of relative
count frequencies in the two groups being compared.}
\item{expression.levels}{a matrix of estimated gene
expression levels as indicated by mean relative read
frequencies. It has three columns \code{grp1}, \code{grp2},
\code{pooled} corresponding to the two treatment groups and
the pooled mean.} \item{log.fc}{base 2 log fold change in
mean relative frequency between two groups.}
\item{p.values}{p-values of the exact NB test applied to
each gene (row).} \item{q.values}{q-values (estimated FDR)
corresponding to the p-values.}
}
\description{
\code{exact.nb.test} performs the Robinson and Smyth exact
negative binomial (NB) test for differential gene
expression on each gene and summarizes the results using
p-values and q-values (FDR).
}
\details{
The negative binomial (NB) distribution offers a more
realistic model for RNA-Seq count variability and still
permits an exact (non-asymptotic) test for comparing
expression levels in two groups.

For each gene, let \eqn{S_1}, \eqn{S_2} be the sums of gene
counts from all biological replicates in each group. The
exact NB test is based on the conditional distribution of
\eqn{S_1|S_1+S_2}: a value of \eqn{S_1} that is too big or
too small, relative to the sum \eqn{S_1+S_2}, indicates
evidence for differential gene expression.  When the
effective library sizes are the same in all replicates and
the dispersion parameters are known, we can determine the
probability functions of \eqn{S_1}, \eqn{S_2} explicitly.
The exact p-value is computed as the total conditional
probability of all possible values of \eqn{(S_1, S_2)} that
have the same sum as but are more extreme than the observed
values of \eqn{(S_1, S_2)}.

Note that we assume that the NB dispersion parameters for
the two groups are the same and library sizes (column
totals of the count matrix) are the same.
}
\note{
Before calling \code{\link{exact.nb.test}}, the user should
call \code{\link{estimate.norm.factors}} to estimate
normalization factors, call \code{\link{prepare.nbp}} to
adjust library sizes, and call \code{\link{estimate.disp}}
to fit a dispersion model. The exact NB test will be
performed using \code{pseudo.counts} in the list
\code{obj}, which are normalized and adjusted to have the
same effective library sizes (column sums of the count
matrix, multiplied by normalization factors).

Users not interested in fine tuning the underlying
statistical model should use \code{\link{nbp.test}}
instead. The all-in-one function \code{\link{nbp.test}}
uses sensible approaches to normalize the counts, estimate
the NBP model parameters and test for differential gene
expression.

A test will be performed on a row (a gene) only when the
total row count is nonzero, otherwise an NA value will be
assigned to the corresponding p-value and q-value.
}
\examples{
## Load Arabidopsis data
data(arab);

## Specify treatment groups
## grp.ids = c(1, 1, 1, 2, 2, 2); # Numbers or strings are both OK
grp.ids = rep(c("mock", "hrcc"), each=3);

## Estimate normalization factors
norm.factors = estimate.norm.factors(arab);
print(norm.factors);

## Prepare an NBP object, adjust the library sizes by thinning the
## counts. For demonstration purpose, only use the first 100 rows of
## the arab data.
set.seed(999);
obj = prepare.nbp(arab[1:100,], grp.ids, lib.size=colSums(arab), norm.factors=norm.factors);
print(obj);

## Fit a dispersion model (NBQ by default)
obj = estimate.disp(obj);
plot(obj);

## Perform exact NB test
## grp1 = 1;
## grp2 = 2;
grp1 = "mock";
grp2 = "hrcc";

obj = exact.nb.test(obj, grp1, grp2);

## Print and plot results
print(obj);
par(mfrow=c(3,2));
plot(obj);
}
\seealso{
\code{\link{nbp.test}}.
}

