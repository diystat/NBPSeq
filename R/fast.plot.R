binxy = function(x, y, xlim=range(x), ylim=range(y), nbins, drop=TRUE) {

  eps = 1e-2/nbins;
  xlim[2] = xlim[2] + diff(xlim)*eps;
  ylim[2] = ylim[2] + diff(ylim)*eps;

  ## Map x to integers between 0 and (nbins-1);
  xf = nbins/diff(xlim);
  yf = nbins/diff(ylim);
  
  xi = floor((x - xlim[1])*xf);
  yi = floor((y - ylim[1])*yf);
  freqs = tabulate(xi + yi*nbins+1, nbins*nbins);

  x = rep(((1:nbins) - 0.5)/xf + xlim[1], nbins);
  y = rep(((1:nbins) - 0.5)/yf + ylim[1], each=nbins);

  ## 
  id = if (drop) {
    freqs > 0
  } else {
    TRUE
  }

  data.frame(x = x[id], y=y[id], freqs= freqs[id]);
}


example.table = function() {
  rm(list=ls());
  source('fast.plot.R');

  set.seed(999);
  x = rnorm(20000);
  y = rnorm(20000, sd=2);
  nbins = resolution = 50;
  h = hist2d(x, y, nbins=nbins);

  ## image(h);
  ## h$z[h$z==0] = NA;

  nblack = 95;
  nred = 5;

  clip = 20;
  col = gray( (200:0)/256);

  image(h, col = col, zlim=c(0.1, clip));

  if (clip < max(h$z)) {
    col = rgb((0:256)/256, 0, 0);
    image(h, col = col, zlim=c(clip, max(h$z)), add=TRUE);
  }

  he = h;
  he$x = exp(h$x);
  he$y = exp(h$y);
  image(he, col = col, zlim=c(0.1, clip), log="xy");


  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((33:256)/256);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }


  t4 = binxy(x, y, xlim = c(-5, 5), ylim=c(-5, 5), nbins=resolution);
  plot(t4$x, t4$y);

  microbenchmark(
    image(h),
    plot(t4$x, t4$y)
    );



  points(x[1:10], y[1:10]);

  xlim = range(x);
  ylim = range(y);
  resolution = 50;
  eps = 1e-2/resolution;
  xlim[2] = xlim[2] + diff(xlim)*eps;
  ylim[2] = ylim[2] + diff(ylim)*eps;

  xi = floor((x - xlim[1])/diff(xlim)*resolution);
  yi = floor((y - ylim[1])/diff(ylim)*resolution);

  ## xi = as.integer((x - xlim[1])/diff(xlim)*resolution);
  ## yi = as.integer((y - ylim[1])/diff(ylim)*resolution);

  ## t1 = table(factor(xi, levels=0:99), factor(yi, levels=0:99));
  t0 = table(xi, yi);
  t2 = tabulate(xi + yi*resolution+1, resolution*resolution);
  t3 = binxy(x, y, nbins=resolution);

  t4 = binxy(x, y, xlim = c(-5, 5), ylim=c(-5, 5), nbins=resolution);
  plot(t4$x, t4$y);

  t4 = binxy(x, y, xlim = c(-5, 5), ylim=c(-5, 5), nbins=resolution, drop=FALSE);
  plot(t4$x, t4$y);

  z = matrix(t4$freqs, resolution, resolution);

  ## .External.graphics(C_image, t4$x, t4$y, t4$freqs, 1);
  ##?.External.graphics

  image(t4$x[1:resolution], t4$y[1:resolution], matrix(r4$z, resolution, resolution));

  ## t4$z = t4$freqs;

  t1 = table(factor(xi, levels=0:(resolution-1)), factor(yi, levels=0:(resolution-1)));
  t2 = tabulate(xi + yi*resolution+1, resolution*resolution);
  dim(t2) = c(resolution, resolution);
 range(t1 - t2);
  
  library(microbenchmark);
  microbenchmark(
  t0 = table(xi, yi),
  t1 = table(factor(xi, levels=0:(resolution-1)), factor(yi, levels=0:(resolution-1))),
  t2 = tabulate(xi + yi*resolution+1, resolution*resolution)
    );

}

example.his2d.old = function() {
  source('fast.plot.R');

  set.seed(999);
  x = rnorm(10);
  y = rnorm(10);
  debug(smart.plot);
  smart.plot(x, y, resolution=100, pch=19, plot=FALSE, log="x");

  smart.plot(x, y, resolution=100, pch=19);

  x = rnorm(20000);
  y = rnorm(20000);

  smart.plot(x, y, resolution=40, pch=19, clip=16);

  obj = smart.plot(x, y, resolution=100, pch=19, plot=FALSE);
  sum(obj$freq);

  smart.plot(x, y, resolution=100, pch=19, log="x");
  smart.plot(x, y, resolution=20, pch=19, log="x", clip=16);

  obj = smart.plot(x, y, resolution=20, pch=19, plot=FALSE, log="x", clip=16);
  sum(x>0);
  sum(obj$freq);
  
  ## debug(smart.plot);
  ## debug(plot.default);
  ## undebug(plot.default);
  par(mfrow=c(2,3));
  plot(x, y, pch=19);
  ## Plot with translucent color
  plot(x, y, col=rgb(0, 0, 1, 0.1), pch=19);
  smart.plot(x, y, resolution=100, pch=19);


  smart.plot(x, y, pch=19);
  smart.plot(x, y, resolution = 100, pch="+");
  smart.plot(x, y, resolution = 100, pch="+", col=1);

  smart.plot(x, y, resolution=100, clip=5, pch=19);
  smart.plot(x, y, resolution= 100, pch=19, log="xy");

  obj = smart.plot(x, y, resolution= 100, pch=19);

  obj = smart.plot(x, y, resolution= 100, pch=19, plot=FALSE);


}
