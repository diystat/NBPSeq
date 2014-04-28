\name{make.disp}
\alias{make.disp}
\title{(private) Specify a dispersion model}
\usage{
make.disp(nb.data, x, model, predictor, subset = filter.mu.pre(nb.data, x),
  ...)
}
\arguments{
  \item{nb}{NB data, output from
  \code{\link{prepare.nb.data}}}

  \item{x}{a matrix, design matrix (specifying the
  treatment structure).}

  \item{model}{a string giving the name of the disperion
  model, can be one of "NB2", "NBP", "NBQ", "NBS" or "step"
  (not case sensitive).}

  \item{predictor}{a string giving the name of the
  predictor to use in the dispersion model, can be one of
  "pi" and "mu", or "rs". \code{"pi"}, preliminarily
  estimated mean relative frequencies; \code{"mu"},
  preliminarily estimated mean frequencies; \code{"rs"},
  row sums.}

  \item{subset}{a list of logical,}

  \item{...}{additional parameter to \code{disp.fun.*}}
}
\value{
a list, output from the call to the funtion
\code{disp.fun.<model>}.
}
\description{
Specify a dispersion model. The parameters of the specified
model are to be estimated from the data using the function
\code{optim.disp.apl} or \code{optim.disp.pl}.
}
\details{
This functions calls \code{disp.fun.<model>} to specify a
dispersion model (a list), using output from a call to
\code{disp.predictor.<predictor>} as argument list, where
\code{<model>} is \code{model} from the input in lower case
(one of "nb2", "nbp", "nbq", "nbs" or "step") and
\code{<predictor>} is \code{predictor} from the input (one
of "pi", "mu", or "rs")
}

