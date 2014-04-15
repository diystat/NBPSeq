\name{test.coefficient}
\alias{test.coefficient}
\title{Large-sample Test for a Regression Coefficient in an
Negative Binomial Regression Model}
\usage{
test.coefficient(nb, dispersion, x, beta0, tests = c("HOA", "LR", "Wald"),
  alternative = "two.sided", subset = 1:m, print.level = 1)
}
\arguments{
  \item{nb}{an NB data object, output from
  \code{\link{prepare.nb.data}}.}

  \item{dispersion}{dispersion estimates, output from
  \code{\link{estimate.disp}}.}

  \item{x}{an \eqn{n} by \eqn{p} design matrix describing
  the treatment structure}

  \item{beta0}{a \eqn{p}-vector specifying the null
  hypothesis. Non-NA components specify the parameters to
  test and their null values. (Currently, only
  one-dimensional test is implemented, so only one non-NA
  component is allowed).}

  \item{tests}{a character string vector specifying the
  tests to be performed, can be any subset of \code{"HOA"}
  (higher-order asymptotic test), \code{"LR"} (likelihood
  ratio test), and \code{"Wald"} (Wald test).}

  \item{alternative}{a character string specifying the
  alternative hypothesis, must be one of \code{"two.sided"}
  (default), \code{"greater"} or \code{"less"}.}

  \item{subset}{an index vector specifying on which rows
  should the tests be performed}

  \item{print.level}{a number controlling the amount of
  messages printed: 0 for suppressing all messages, 1
  (default) for basic progress messages, and 2 to 5 for
  increasingly more detailed message.}
}
\value{
a list containing the following components:
\item{beta.hat}{an \eqn{m} by \eqn{p} matrix of regression
coefficient under the full model} \item{mu.hat}{an \eqn{m}
by \eqn{n} matrix of fitted mean frequencies under the full
model} \item{beta.tilde}{an \eqn{m} by \eqn{p} matrix of
regression coefficient under the null model}
\item{mu.tilde}{an \eqn{m} by \eqn{n} matrix of fitted mean
frequencies under the null model.} \item{HOA, LR,
Wald}{each is a list of two \eqn{m}-vectors,
\code{p.values} and \code{q.values}, giving p-values and
q-values of the corresponding tests when that test is
included in \code{tests}.}
}
\description{
\code{test.coefficient} performs large-sample tests
(higher-order asymptotic test, likelihood ratio test,
and/or Wald test) for testing regression coefficients in an
NB regression model.
}
\details{
\code{test.coefficient} performs large-sample tests for a
one-dimensional (\eqn{q=1}) component \eqn{\psi} of the
\eqn{p}-dimensional regression coefficient \eqn{\beta}. The
hypothesized value \eqn{\psi_0} of \eqn{\psi} is specified
by the non-NA component of the vector \code{beta0} in the
input.

The likelihood ratio statistic, \deqn{ \lambda = 2
(l(\hat\beta) - l(\tilde\beta)),} converges in distribution
to a chi-square distribution with \eqn{1} degree of
freedom.  The signed square root of the likelihood ratio
statistic \eqn{\lambda}, also called the directed deviance,
\deqn{r = sign (\hat\psi - \psi_0) \sqrt \lambda} converges
to a standard normal distribution.

For testing a one-dimensional parameter of interest,
Barndorff-Nielsen (1986, 1991) showed that a modified
directed \deqn{ r^* = r - \frac{1}{r} \log(z)} is, in wide
generality, asymptotically standard normally distributed to
a higher order of accuracy than the directed deviance
\eqn{r} itself, where \eqn{z} is an adjustment term. Tests
based on high-order asymptotic adjustment to the likelihood
ratio statistic, such as \eqn{r^*} or its approximation,
are referred to as higher-order asymptotic (HOA) tests.
They generally have better accuracy than corresponding
unadjusted likelihood ratio tests, especially in situations
where the sample size is small and/or when the number of
nuisance parameters (\eqn{p-q}) is large. The
implementation here is based on Skovgaard (2001). See Di et
al. 2013 for more details.
}
\examples{
## Load Arabidopsis data
data(arab);

## Estimate normalization factors (we want to use the entire data set)
norm.factors = estimate.norm.factors(arab);

## Prepare the data
## For demonstration purpose, only the first 50 rows are used
nb.data = prepare.nb.data(arab[1:50,], lib.sizes = colSums(arab), norm.factors = norm.factors);

## For real analysis, we will use the entire data set, and can omit lib.sizes parameter)
## nb.data = prepare.nb.data(arab, norm.factors = norm.factors);

print(nb.data);
plot(nb.data);

## Specify the model matrix (experimental design)
grp.ids = as.factor(c(1, 1, 1, 2, 2, 2));
x = model.matrix(~grp.ids);

## Estimate dispersion model
dispersion = estimate.dispersion(nb.data, x);

print(dispersion);
plot(dispersion);

## Specify the null hypothesis
## The null hypothesis is beta[2]=0 (beta[2] is the log fold change).
beta0 = c(NA, 0);

## Test regression coefficient
res = test.coefficient(nb.data, dispersion, x, beta0);

## The result contains the data, the dispersion estimates and the test results
print(str(res));

## Show HOA test results for top ten most differentially expressed genes
top = order(res$HOA$p.values)[1:10];
print(cbind(nb.data$counts[top,], res$HOA[top,]));

## Plot log fold change versus the fitted mean of sample 1 (analagous to an MA-plot).
plot(res$mu.tilde[,1], res$beta.hat[,2]/log(2), log="x",
     xlab="Fitted mean of sample 1 under the null",
     ylab="Log (base 2) fold change");

## Highlight top DE genes
points(res$mu.tilde[top,1], res$beta.hat[top,2]/log(2), col="magenta");



}
\references{
Barndorff-Nielsen, O. (1986): "Infereni on full or partial
parameters based on the standardized signed log likelihood
ratio," Biometrika, 73, 307-322

Barndorff-Nielsen, O. (1991): "Modified signed log
likelihood ratio," Biometrika, 78, 557-563.

Skovgaard, I. (2001): "Likelihood asymptotics,"
Scandinavian Journal of Statistics, 28, 3-32.

Di Y, Schafer DW, Emerson SC, Chang JH (2013): "Higher
order asymptotics for negative binomial regression
inferences from RNA-sequencing data". Stat Appl Genet Mol
Biol, 12(1), 49-70.
}

