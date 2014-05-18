##' Commpute a 2-d histogram of the given data values. (not implemented yet: If \code{plot == TRUE}, plot the resulting histogram.)
##'
##' This funciton divides the \code{xlim x ylim} region into \code{nbins x nbins}
##' equal-sized cells and count the number of (x,y) points in each cell.
##' @note
##'
##' Only points inside the region defined by \code{xlim x ylim}
##' (inclusive) will be counted. For each cell, the lower boundaries
##' are closed and upper boundaries are open.  A small number will be
##' added to the upper limits in xlim and ylim so that no points will
##' be on the region's upper boundaries.
##'
##' @title 2-d Histogram
##' @param x a vector
##' @param y a vector of the same length as x
##' @param xlim a vector of length 2, the range of x values 
##' @param ylim a vector of length 2, the range of y values
##' @param nbins a single number giving the number of bins (the same for both x- and y- axes).
##' @return a list
##' \item{x}{a vector of length nbins, the midpoints of each bin on the x-axis.}
##' \item{y}{a vector of length nbins, the midpoints of each bin on the y-axis.}
##' \item{z}{a \code{nbins} by \code{nbins} matrix of of counts. For each cell, the number of (x, y) inside.}
##' The list can be passed to \code{image()} directly for potting.
hist2d= function(x, y, xlim=range(x), ylim=range(y), nbins) {

  id = (x >= xlim[1]) & (x <= xlim[2]) & (y >=ylim[1]) & (y <= ylim[2]);
  x = x[id];
  y = y[id];

  eps = 1e-5/nbins;
  xlim[2] = xlim[2] + diff(xlim)*eps;
  ylim[2] = ylim[2] + diff(ylim)*eps;

  ## Map x to integers between 0 and (nbins-1);
  xf = nbins/diff(xlim);
  yf = nbins/diff(ylim);
  
  xi = floor((x - xlim[1])*xf);
  yi = floor((y - ylim[1])*yf);
  freqs = tabulate(xi + yi*nbins+1, nbins*nbins);

  list(x = ((1:nbins) - 0.5)/xf + xlim[1],
       y = ((1:nbins) - 0.5)/yf + ylim[1],
       z = matrix(freqs, nbins, nbins));
}

example.hist2d = function() {
  source('smart.plot.R');
  x = 1:10
  y = (1:10) * 10; 

  par(mfrow=c(2,2));
  ## The points may not be in the center of each cell
  h = hist2d(x, y, nbins=10);
  image(h);
  points(x, y);

  ## If that matters, choose  xlim and ylim more carefully
  h = hist2d(x, y, xlim=c(0.5, 10.5), ylim=c(5, 105), nbins=10);
  image(h);
  points(x, y);

  ## For higher resolution, this shouldn't be a problem
  h = hist2d(x, y, nbins=50);
  image(h);
  points(x, y);

  ## Note in plot the x and y axes are extended 4% beyond the xlim and ylim by default.
  plot(x, y);

  ## Color choice when using image
  n = 1e5;
  x = rnorm(n);
  y = rnorm(n);
  h = hist2d(x,y, nbins=20);
  h$z[h$z==0] = NA;
  print(h);

  par(mfrow=c(2,2));
  image(h, col=gray((256:0)/256));
  image(h, col=gray((0:256)/256));
  image(h, col=gray((224:0)/256));

  ## Adding contour lines
  par(mfrow=c(1,1)); 
  image(h, col=gray((224:0)/256));
  contour(h, nlevels=5, col="magenta", add=TRUE);

  ## We can not use Inf as zlim
  image(h, zlim=c(0, Inf));
  image(h, zlim=c(0, 100));
  image(h, zlim=c(0, 1e5));

  ## 
  n = 1e4;
  x = rnorm(n);
  y = rexp(n);
  h = hist2d(x,y, nbins=20);
  h$z = log(h$z);
  print(h);
  image(h);
  points(x, y, pch=4, col="cyan");

  system.time({h = hist2d(x,y, nbins=20)});

  ## Timing
  library(microbenchmark);
  n = 1e4;
  x = rnorm(n);
  y = rexp(n);
  h = hist2d(x,y, nbins=20);
  microbenchmark(
    h = hist2d(x,y, nbins=20), 
    i = image(h),
    c = contour(h)
    );

}



##' An alternative to plot.default() for plotting a large number of
##' densely distributed points.  This function can produce a visually
##' almost identical plot using only a subset of the points.  This is
##' particular useful for reducing output file size when plots are
##' written to eps files.
##'
##'   Writing plots with a large number of points to eps files can
##'   result in big files and lead to very slow rendering time.
##'
##'   Usually for a large number of points, a lot of them will
##'   overlap with each other. Plotting only a subset of selected
##'   non-overlapping points can give visually almost identical
##'   plots. Further more, the plots can be enhanced if using gray
##'   levels (the default setting) that are proportional to the
##'   number points overlapping with each plotted point.
##' 
##'   This function scans the points sequentially. For each unmarked
##'   point that will be plotted, all points that overlap with it
##'   will be marked and not to plotted, and the number of
##'   overlapping points will be recorded. This is essentially
##'   producing a 2d histogram. The freqs of the points will be
##'   converted to gray levels, darker colors correspond to higher
##'   freqs.
##'
##' @title (private) An alternative to plot.default() for plotting a large number of
##' densely distributed points.
##'
##' @param x  x
##' @param y  y
##' @param xlim xlim
##' @param ylim ylim
##' @param xlab x label
##' @param ylab y label
##' @param log log
##' @param resolution a number, determines the distance below which
##' points will be considered as overlapping.
##' @param plot logical, whether
##' @param col color
##' @param clip clip
##' @param color.clipped color of clipped points
##' @param ... other arguments are the same as in plot.default().
##'
##' @return (if plot=FALSE) a list 
##'
##'   \item{x, y}{the x, y-coordinates of the subset of representative points}
##' 
##'   \item{id}{the indicies of these points in the original data set}
##'  
##'   \item{freqs}{the numbers of points that overlap with each representative point}
##'
##'   \item{col}{colors determined by the freqs}
smart.plot.new = function(x, y=NULL, xlim=NULL, ylim=NULL,
  xlab=NULL, ylab=NULL,
  log="", resolution=50,
  col = gray( (224:0)/256),
  clip=NULL,
  col.clipped = rgb(log2(1:256)/log2(256), 0, 0),
  ## col.clipped = rgb(c(0:256)/256, 0, 0),
  ...) {

  ## These lines are copied from plot.default.
  xlabel = if (!missing(x)) deparse(substitute(x));
  ylabel = if (!missing(y)) deparse(substitute(y));
  xy = xy.coords(x, y, xlabel, ylabel, log);
  xlab = if (is.null(xlab)) xy$xlab else xlab;
  ylab = if (is.null(ylab)) xy$ylab else ylab;
  xlim = if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim;
  ylim = if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim;
  
  n = length(xy$x);

  id = is.finite(xy$x) & is.finite(xy$y);
  ## id[xy$x < xlim[1] | xy$x > xlim[2] | xy$y < ylim[1] | xy$y > ylim[2]]=FALSE;

  logxy = strsplit(log, NULL)[[1L]];
  logx = 'x' %in% logxy;
  logy = 'y' %in% logxy;

  xx = if (logx) log(xy$x[id]) else xy$x[id];
  yy = if (logy) log(xy$y[id]) else xy$y[id];

  xxlim = if (logx) log(xlim) else xlim;
  yylim = if (logy) log(ylim) else ylim;

  h = hist2d(xx, yy, xxlim, yylim, resolution);

  ## When calling image later, we would like to pass the log
  ## parameter. For that to work, we need transform h$x to exp(h$x) if
  ## logx is true.
  if (logx) h$x = exp(h$x);
  if (logy) h$y = exp(h$y);

  ## We want to use two sets of colors, one for frequencies 1:clip;
  ## another for frequencies greater than clip

  ## col = gray( (200:0)/256);
  ## col.clipped = rgb((0:256)/256, 0, 0);

  if (is.null(clip)) {
    image(h, col = col, zlim=c(0.1, max(h$z)), log=log, xlab=xlab, ylab=ylab, ...);
  } else {
    if (clip < max(h$z)) {
      image(h, col = col, zlim=c(0.1, clip), log=log, xlab=xlab, ylab=ylab, ...);
      image(h, col = col.clipped, zlim=c(clip, max(h$z)), add=TRUE, log=log);
    } else {
      image(h, col = col, zlim=c(0.1, clip), log=log, xlab=xlab, ylab=ylab, ...);
    }
  }

  invisible(h)
}

example.smart.plot = function() {
  source('smart.plot.R');

  image(h);
  points(x, y);
  

  set.seed(999);
  n = 10000;
  x = rlnorm(n);
  y = rnorm(n);

  par(mfrow=c(2,2));
  h = smart.plot(x, y, resolution=10);
  h = smart.plot(x, y, resolution=10, clip=20);
  ## debug(smart.plot);
  h = smart.plot(x, y, resolution=10, log="x");
  h = smart.plot(x, y, resolution=10, log="x", clip=20);

  par(mfrow=c(2,2));
  smart.plot(x, y);
  smart.plot(x, y, clip=20);

  ## debug(smart.plot);
  h = smart.plot(x, y, log="x");
  h = smart.plot(x, y, log="x", clip=20);

  library(NBPSeq);
  par(mfrow=c(1,1));
  data(arab);
  h = smart.plot(arab[,1:2], log="xy");
  ## contour(h, add=TRUE, col="cyan", nlevels=5);
  ## h
  
  h = smart.plot(arab[,1:2], log="xy", clip=2000);

  h = smart.plot(arab[,1:2], log="xy", clip=128, resolution=20,
    main="Arabidopsis", xlab="Sample 1", ylab="Sample 2");
  abline(0, 1, col="cyan")
  contour(h, add=TRUE, col="cyan");

}

##' An alternative to plot.default() for plotting a large number of
##' densely distributed points.  This function can produce a visually
##' almost identical plot using only a subset of the points.  This is
##' particular useful for reducing output file size when plots are
##' written to eps files.
##'
##'   Writing plots with a large number of points to eps files can
##'   result in big files and lead to very slow rendering time.
##'
##'   Usually for a large number of points, a lot of them will
##'   overlap with each other. Plotting only a subset of selected
##'   non-overlapping points can give visually almost identical
##'   plots. Further more, the plots can be enhanced if using gray
##'   levels (the default setting) that are proportional to the
##'   number points overlapping with each plotted point.
##' 
##'   This function scans the points sequentially. For each unmarked
##'   point that will be plotted, all points that overlap with it
##'   will be marked and not to plotted, and the number of
##'   overlapping points will be recorded. This is essentially
##'   producing a 2d histogram. The freqs of the points will be
##'   converted to gray levels, darker colors correspond to higher
##'   freqs.
##'
##' @title (private) An alternative to plot.default() for plotting a large number of
##' densely distributed points.
##'
##' @param x  x
##' @param y  y
##' @param xlim xlim
##' @param ylim ylim
##' @param xlab x label
##' @param ylab y label
##' @param log log
##' @param resolution a number, determines the distance below which
##' points will be considered as overlapping.
##' @param plot logical, whether
##' @param col color
##' @param clip clip
##' @param color.clipped color of clipped points
##' @param ... other arguments are the same as in plot.default().
##'
##'
##' @return (if plot=FALSE) a list 
##'
##'   \item{x, y}{the x, y-coordinates of the subset of representative points}
##' 
##'   \item{id}{the indicies of these points in the original data set}
##'  
##'   \item{freqs}{the numbers of points that overlap with each representative point}
##'
##'   \item{col}{colors determined by the freqs}
smart.plot.old = function(x, y=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
  log="", resolution=100, plot=TRUE, col=NULL, clip=Inf, color.clipped=TRUE, ...) {


  ## These lines are copied from plot.default.
  xlabel = if (!missing(x)) deparse(substitute(x));
  ylabel = if (!missing(y)) deparse(substitute(y));
  xy = xy.coords(x, y, xlabel, ylabel, log);
  xlab = if (is.null(xlab)) xy$xlab else xlab;
  ylab = if (is.null(ylab)) xy$ylab else ylab;
  xlim = if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim;
  ylim = if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim;
  
  x = xy$x;
  y = xy$y;
  n = length(x);
  id = is.finite(xy$x) & is.finite(xy$y);
  id[id & (x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])]=FALSE;

  logxy = strsplit(log, NULL)[[1L]];
  if ("x" %in% logxy) {
    x[id] = log(x[id]);
    epsx = diff(log(xlim)) / resolution;
  } else
    epsx = diff(xlim) / resolution;

  if ("y" %in% logxy) {
    y[id] = log(y[id]);
    epsy = diff(log(ylim)) / resolution; 
  } else
    epsy = diff(ylim) / resolution; 

  counts = rep(0, n);
  i = 1;

  ## Scan the points and select non-overlapping ones
  while (i < n) {
    if (id[i]) {
      ids = ((i+1):n)[id[(i+1):n]];
      overlap = ids[abs(x[ids] - x[i]) < epsx & abs(y[ids] - y[i]) < epsy];
      id[overlap] = FALSE;
      counts[i] = length(overlap) + 1;
    }
    i = i + 1;
  }

  ## Sort the data so that points representing more points will be plotted later.
  id = (1:n)[id];
  id = id[order(counts[id])];
  counts = counts[id];

  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((1 - counts.clipped / max(counts.clipped)) * 0.8);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }

  if (plot) {
    plot(x=xy$x[id], y=xy$y[id], xlim=xlim, ylim=ylim, log=log, xlab=xlab, ylab=ylab, col=col, ...);
    invisible();
  }

  invisible(list(x=xy$x[id], y=xy$y[id], freqs=counts, col = col, id=id));
}

smart.plot = smart.plot.new;

##' See description of \code{\link{smart.plot}} for more details.
##'
##' @title (private) An alternative to point.default() for plotting a large number of
##' densely distributed points.
##'
##' @param x x
##' @param y y
##' @param resolution a number, determines the distance below which
##' points will be considered as overlapping.
##' @param col color
##' @param clip clip
##' @param color.clipped color of clipped points
##' @param ... other arguments are the same as in plot.default().
##'
##' @return NULL
smart.points = function(x, y=NULL, 
  resolution=50, col=NULL, clip=Inf, color.clipped=TRUE, ...) {

  xy = xy.coords(x, y);

  xlim = range(xy$x[is.finite(xy$x)]);
  ylim = range(xy$y[is.finite(xy$y)]);

  x = xy$x;
  y = xy$y;
  n = length(x);
  id = is.finite(xy$x) & is.finite(xy$y);
  id[id & (x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])]=FALSE;

  log = "";
  
  logxy = strsplit(log, NULL)[[1L]];
  if ("x" %in% logxy) {
    x[id] = log(x[id]);
    epsx = diff(log(xlim)) / resolution;
  } else
    epsx = diff(xlim) / resolution;

  if ("y" %in% logxy) {
    y[id] = log(y[id]);
    epsy = diff(log(ylim)) / resolution; 
  } else
    epsy = diff(ylim) / resolution; 

  counts = rep(0, n);
  i = 1;

  ## Scan the points and select non-overlapping ones
  while (i < n) {
    if (id[i]) {
      ids = ((i+1):n)[id[(i+1):n]];
      overlap = ids[abs(x[ids] - x[i]) < epsx & abs(y[ids] - y[i]) < epsy];
      id[overlap] = FALSE;
      counts[i] = length(overlap) + 1;
    }
    i = i + 1;
  }

  ## Sort the data so that points representing more points will be plotted later.
  id = (1:n)[id];
  id = id[order(counts[id])];
  counts = counts[id];

  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((1 - counts.clipped / max(counts.clipped)) * 0.8);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }

  points(x=xy$x[id], y=xy$y[id], col=col, ...);
  invisible();
}

example.smart.plot = function() {
  x = rnorm(20000);
  y = rnorm(20000);

  ## debug(smart.plot);
  ## debug(plot.default);
  ## undebug(plot.default);
  par(mfrow=c(2,3));
  plot(x, y, pch=19);
  ## Plot with translucent color
  plot(x, y, col=rgb(0, 0, 1, 0.1), pch=19);
  smart.plot(x, y, resolution=100);

  smart.plot(x, y, pch=19);
  smart.plot(x, y, resolution = 100, pch="+");
  smart.plot(x, y, resolution = 100, pch="+", col=1);

  smart.plot(x, y, resolution=100, clip=5, pch=19);
  smart.plot(x, y, resolution= 100, pch=19, log="xy");

  obj = smart.plot(x, y, resolution= 100, pch=19);

##  obj = smart.plot(x, y, resolution= 100, pch=19, plot=FALSE);


}
