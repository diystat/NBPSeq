\name{estimate.disp}
\alias{estimate.disp}
\title{Fit a parametric disperison model to thinned counts}
\usage{
estimate.disp(obj, model = "NBQ", print.level = 1, ...)
}
\arguments{
  \item{obj}{output from \code{\link{prepare.nbp}}.}

  \item{model}{a string, one of "NBQ" (default), "NBP" or
  "NB2".}

  \item{print.level}{a number, controls the amount of
  messages printed: 0 for suppressing all messages, 1 for
  basic progress messages, larger values for more detailed
  messages.}

  \item{\dots}{additional parameters controlling the
  estimation of the parameters.}
}
\value{
The list \code{obj} from the input with some added
components summarizing the fitted dispersion model.  Users
can print and plot the output to see brief summaries of the
fitted dispersion model. The output is otherwise not
intended for use by end users directly.
}
\description{
Fit a parametric dispersion model to RNA-Seq counts data
prepared by \code{\link{prepare.nbp}}. The model parameters
are estimated from the pseudo counts: thinned/down-sampled
counts that have the same effective library size.
}
\details{
For each individual gene \eqn{i}, a negative binomial (NB)
distribution uses a dispersion parameter \eqn{\phi_i} to
capture the extra-Poisson variation between biological
replicates: the NB model imposes a mean-variance
relationship \eqn{\sigma_i^2 = \mu_i + \phi_i \mu_i^2}.  In
many RNA-Seq data sets, the dispersion parameter
\eqn{\phi_i} tends to vary with the mean \eqn{\mu_i}. We
proposed to capture the dispersion-mean dependence using
parametric models.

With this function, \code{estimate.disp}, users can choose
from three parametric models: NB2, NBP and NBQ (default).

Under the NB2 model, the dispersion parameter is a constant
and does not vary with the mean expression levels.

Under the NBP model, the log dispersion is modeled as a
linear function of preliminarily estimated log mean
relative frequencies (\code{pi.pre}):

log(phi) = par[1] + par[2] * log(pi.pre/pi.offset),

Under the NBQ model, the log dispersion is modeled as a
quadratic function of preliminarily estimated log mean
relative frequencies (\code{pi.pre}):

log(phi) = par[1] + par[2] * log(pi.pre/pi.offset) + par[3]
* (log(pi.pre/pi.offset))^2;

The NBQ model is more flexible than the NBP and NB2 models,
and is the current default option.

In the NBP and NBQ models, \code{pi.offset} is fixed to be
1e-4, so par[1] corresponds to the dispersion level when
the relative mean frequency is 100 reads per million (RPM).

The dispersion parameters are estimated from the pseudo
counts (counts adjusted to have the same effective library
sizes).  The parameters are estimated by maximizing the log
conditional likelihood of the model parameters given the
row sums. The log conditional likelihood is computed for
each gene in each treatment group and then summed over
genes and treatment groups.
}
\note{
Users should call \code{\link{prepare.nbp}} before calling
this function. The function \code{\link{prepare.nbp}} will
normalize the counts and adjust the counts so that the
effective library sizes are approximately the same
(computing the conditional likelihood requires the library
sizes to be the same).
}
\examples{
## See the example for nb.exact.test
}
\references{
Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP
Negative Binomial Model for Assessing Differential Gene
Expression from RNA-Seq", Statistical Applications in
Genetics and Molecular Biology, 10 (1).
}
\seealso{
\code{\link{nbp.test}}, \code{\link{exact.nb.test}}
}

