\name{smart.plot}
\alias{smart.plot}
\title{(private) An alternative to plot.default() for plotting a large number of
densely distributed points.}
\usage{
smart.plot(x, y = NULL, xlim = NULL, ylim = NULL, xlab = NULL,
  ylab = NULL, log = "", resolution = 100, plot = TRUE, col = NULL,
  clip = Inf, color.clipped = TRUE, ...)
}
\arguments{
  \item{x}{x}

  \item{y}{y}

  \item{xlim}{xlim}

  \item{ylim}{ylim}

  \item{xlab}{x label}

  \item{ylab}{y label}

  \item{log}{log}

  \item{resolution}{a number, determines the distance below
  which points will be considered as overlapping.}

  \item{plot}{logical, whether}

  \item{col}{color}

  \item{clip}{clip}

  \item{color.clipped}{color of clipped points}

  \item{...}{other arguments are the same as in
  plot.default().}
}
\value{
(if plot=FALSE) a list

\item{x, y}{the x, y-coordinates of the subset of
representative points}

\item{id}{the indicies of these points in the original data
set}

\item{freqs}{the numbers of points that overlap with each
representative point}

\item{col}{colors determined by the freqs}
}
\description{
An alternative to plot.default() for plotting a large
number of densely distributed points.  This function can
produce a visually almost identical plot using only a
subset of the points.  This is particular useful for
reducing output file size when plots are written to eps
files.
}
\details{
Writing plots with a large number of points to eps files
can result in big files and lead to very slow rendering
time.

Usually for a large number of points, a lot of them will
overlap with each other. Plotting only a subset of selected
non-overlapping points can give visually almost identical
plots. Further more, the plots can be enhanced if using
gray levels (the default setting) that are proportional to
the number points overlapping with each plotted point.

This function scans the points sequentially. For each
unmarked point that will be plotted, all points that
overlap with it will be marked and not to plotted, and the
number of overlapping points will be recorded. This is
essentially producing a 2d histogram. The freqs of the
points will be converted to gray levels, darker colors
correspond to higher freqs.
}

