##' @title Boxplot and scatterplot matrix of relative frequencies (after normalization)
##' @export
##'
##' @param x output from \code{\link{prepare.nb.data}}
##' @param resolution 
##' @param hlim a single number controls the height of the bars in the
##' @param clip 
##' @param eps a small positive number added to rpm
##' @param ... currently not used
##' histograms
##' @return  NULL
plot.nb.data = function(x, resolution=50, hlim =0.25, 
  clip=128, 
  eps = 1e-2,
  ... ) {

  ## 2013-11-22
  ## 2014-05-17, YD

  nb.data = x;

  ## Currently, log is fixed to "xy"
  log = "xy";

  logxy = strsplit(log, NULL)[[1L]];
  logx = 'x' %in% logxy;
  logy = 'y' %in% logxy;

  ## A small positive number

  ## Rows with nonzero reads
  ## id.nonzero = rowSums(nb.data$counts)>0;

  ## boxplot(log(nb.data$rel.freq[id.nonzero,] + eps));
  ## hist(log(nb.data$rel[,1] + eps));

  ## put histograms on the diagonal
  ## used code from the internet by Duncan Murdoch
  panel.hist = function(x, ...) {

    usr = par("usr");
    on.exit(par(usr));
    par(usr = c(usr[1:2], 0, hlim));
    ## print(usr);

    if (logx) {
      h = hist(log(x), plot = FALSE)
      ## print(h);
      breaks = exp(h$breaks);
      yb = 1;
      ## Area corresponding to each bar
      y = exp(diff(h$breaks) * h$density)
      ## points(exp(h$mids), exp(h$density), type="h");
    } else {
      h = hist(x, plot = FALSE)
      ## print(h);
      breaks =h$breaks;
      yb = 0;
      y = diff(h$breaks) * h$density;
      ## y = h$counts/max(h$counts);
    }

    nB = length(breaks)
    ## print(breaks);
    ## print(y);
    rect(breaks[-nB], yb, breaks[-1], y, col="cyan", ...)
  }
  
  pairs(nb.data$rel * 1e6 + 1e-2,
        panel=function(x, ...){
          smart.plot(x, clip=clip, resolution=resolution, add=TRUE, log=log, ...);
          abline(a=0, b=1, col="cyan")
        },
        diag.panel=panel.hist,
        log=log,
        ...
        );

  invisible();
}


