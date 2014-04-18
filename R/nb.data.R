##' @title hello
##' @param x 
##' @param i 
##' @param j 
##' @param ... 
##' @param drop 
##' @return  NULL
`[.nb.data` = function(x, i, j, ..., drop=FALSE) {
  y = list();
  y$counts = x$counts[i,j, ..., drop=drop];
  y$lib.sizes = x$lib.sizes[j];
  y$norm.factors = x$norm.factors[j];
  y$eff.lib.sizes = x$eff.lib.sizes[j];
  y$rel.frequencies = x$rel.frequencies[i,j,..., drop=drop];
  y$tags = x$tags[i,];

  class(y) = class(x);
  y
}

if (FALSE) {
nb[1:10,1:3];
nb[1:10,4:6];
data(arab);
nb = prepare.nb.data(arab);
class(nb) = "nb.data";
str(nb);
}


##' @title Print summary of the nb counts
##' @param x output from \code{\link{prepare.nb.data}}
##' @param ... additional parameters, currently not used
##' @return NULL
print.nb.data = function(x, ...) {

  print(str(x));
  print("Counts:");
  print(head(x$counts));
  cat("...\n");
  print("Lirary sizes (unnormalized):");
  print(x$lib.sizes);
  print("Normalization factors:");
  print(x$norm.factors);
  print("Effective (normalized) library sizes:");
  print(x$eff.lib.sizes);
  print("Reads per Million:");
  print(head(x$rel.freq * 1e6));
  cat("...\n");

  invisible();
}

