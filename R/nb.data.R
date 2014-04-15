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

