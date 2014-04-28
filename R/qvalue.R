qvalue <- function(p=NULL, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, 
  smooth.df = 3, smooth.log.pi0 = FALSE) {
    gui=FALSE;
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#fdr.level: a level at which to control the FDR (optional)
#lambda: the value of the tuning parameter to estimate pi0 (optional)
#pi0.method: either "smoother" or "bootstrap"; the method for automatically
#           choosing tuning parameter in the estimation of pi0, the proportion
#           of true null hypotheses
#robust: an indicator of whether it is desired to make the estimate more robust
#        for small p-values and a direct finite sample estimate of pFDR (optional)
#gui: A flag to indicate to 'qvalue' that it should communicate with the gui.  ## change by Alan
#     Should not be specified on command line.
#smooth.df: degrees of freedom to use in smoother (optional)
#smooth.log.pi0: should smoothing be done on log scale? (optional)
#
#Output
#=============================================================================
#call: gives the function call
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if fdr.level is specified, an indicator of whether the q-value
#    fell below fdr.level (taking all such q-values to be significant controls
#    FDR at level fdr.level)

#Set up communication with GUI, if appropriate
#    print(sys.calls())
#    print(sys.frames())

#    if(gui) {
#        idx <- (1:sys.nframe())[as.character(sys.calls()) == "qvalue.gui()"]
#        gui.env <- sys.frames()[[idx]]
#    }

#This is just some pre-processing
    if(is.null(p))  ## change by Alan
      {qvalue.gui(); return("Launching point-and-click...")}
    if(gui & !interactive())  ## change by Alan
      gui = FALSE

    if(min(p)<0 || max(p)>1) {
      if(gui) ## change by Alan:  check for GUI
        eval(expression(postMsg(paste("ERROR: p-values not in valid range.", "\n"))), parent.frame())
      else
        print("ERROR: p-values not in valid range.")
      return(0)
    }
    if(length(lambda)>1 && length(lambda)<4) {
      if(gui)
        eval(expression(postMsg(paste("ERROR: If length of lambda greater than 1, you need at least 4 values.",
            "\n"))), parent.frame())
      else
        print("ERROR: If length of lambda greater than 1, you need at least 4 values.")
      return(0)
    }
    if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) { ## change by Alan:  check for valid range for lambda
      if(gui)
        eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", "\n"))), parent.frame())
      else
        print("ERROR: Lambda must be within [0, 1).")
      return(0)
    }
    m <- length(p)
#These next few functions are the various ways to estimate pi0
    if(length(lambda)==1) {
        if(lambda<0 || lambda>=1) { ## change by Alan:  check for valid range for lambda
          if(gui)
            eval(expression(postMsg(paste("ERROR: Lambda must be within [0, 1).", "\n"))), parent.frame())
          else
            print("ERROR: Lambda must be within [0, 1).")
          return(0)
        }

        pi0 <- mean(p >= lambda)/(1-lambda)
        pi0 <- min(pi0,1)
    }
    else {
        pi0 <- rep(0,length(lambda))
        for(i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
        }

        if(pi0.method=="smoother") {
            if(smooth.log.pi0)
              pi0 <- log(pi0)

            spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
            pi0 <- predict(spi0,x=max(lambda))$y

            if(smooth.log.pi0)
              pi0 <- exp(pi0)
            pi0 <- min(pi0,1)
        }
        else if(pi0.method=="bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0,length(lambda))
            pi0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (pi0.boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])
            pi0 <- min(pi0,1)
        }
        else {  ## change by Alan: check for valid choice of 'pi0.method' (only necessary on command line)
            print("ERROR: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if(pi0 <= 0) {
      if(gui)
        eval(expression(postMsg(
            paste("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.",
                "\n"))), parent.frame())
      else
        print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
      return(0)
    }
    if(!is.null(fdr.level) && (fdr.level<=0 || fdr.level>1)) {  ## change by Alan:  check for valid fdr.level
      if(gui)
        eval(expression(postMsg(paste("ERROR: 'fdr.level' must be within (0, 1].", "\n"))), parent.frame())
      else
        print("ERROR: 'fdr.level' must be within (0, 1].")
      return(0)
    }
#The estimated q-values calculated here
    u <- order(p)

    # change by Alan
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
 
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
    }

    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
#The results are returned
    if(!is.null(fdr.level)) {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, fdr.level=fdr.level, ## change by Alan
          significant=(qvalue <= fdr.level), lambda=lambda)
    }
    else {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, lambda=lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}

qplot <- function(qobj, rng=c(0.0, 0.1), smooth.df = 3, smooth.log.pi0 = FALSE, ...) { ## change by Alan:  
##  'rng' a vector instead of an upper bound alone
#Input
#=============================================================================
#qobj: a q-value object returned by the qvalue function
#rng: the range of q-values to be plotted (optional)
#smooth.df: degrees of freedom to use in smoother (optional)
#smooth.log.pi0: should smoothing be done on log scale? (optional)
#
#Output
#=============================================================================
#Four plots:
#Upper-left: pi0.hat(lambda) versus lambda with a smoother
#Upper-right: q-values versus p-values
#Lower-left: number of significant tests per each q-value cut-off
#Lower-right: number of expected false positives versus number of significant tests
q2 <- qobj$qval[order(qobj$pval)]
if(min(q2) > rng[2]) {rng <- c(min(q2), quantile(q2, 0.1))} ## change by Alan:  replace 'rng' with vector
p2 <- qobj$pval[order(qobj$pval)]
par(mfrow=c(2,2))
lambda <- qobj$lambda
if(length(lambda)==1) {lambda <- seq(0,max(0.90,lambda),0.05)}
pi0 <- rep(0,length(lambda))
for(i in 1:length(lambda)) {
    pi0[i] <- mean(p2>lambda[i])/(1-lambda[i])
    }
    
if(smooth.log.pi0)
  pi0 <- log(pi0)
spi0 <- smooth.spline(lambda,pi0,df=smooth.df)

if(smooth.log.pi0) {
  pi0 <- exp(pi0)
  spi0$y <- exp(spi0$y)
}

pi00 <- round(qobj$pi0,3)
plot(lambda,pi0,xlab=expression(lambda),ylab=expression(hat(pi)[0](lambda)),pch=".")
mtext(substitute(hat(pi)[0] == that, list(that= pi00)))
lines(spi0)

plot(p2[q2 >= rng[1] & q2 <= rng[2]], q2[q2 >= rng[1] & q2 <= rng[2]], type = "l", xlab = "p-value", ## changes by Alan
  ylab = "q-value")
plot(q2[q2 >= rng[1] & q2 <= rng[2]], (1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), type="l",
  xlab="q-value cut-off", ylab="significant tests")
plot((1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), q2[q2 >= rng[1] & q2 <= rng[2]] *
  (1 + sum(q2 < rng[1])):sum(q2 <= rng[2]), type = "l", xlab = "significant tests",
  ylab = "expected false positives")
par(mfrow=c(1,1))
}

plot.qvalue <- function(x, ...) qplot(x, ...)

qwrite <- function(qobj, filename="my-qvalue-results.txt") {
#Input
#=============================================================================
#qobj: a q-value object returned by the qvalue function
#filename: the name of the file where the results are written
#
#Output
#=============================================================================
#A file sent to "filename" with the following:
#First row: the estimate of the proportion of true negatives, pi0
#Second row: FDR significance level (if specified) ## change by Alan
#Third row and below: the p-values (1st column), the estimated q-values (2nd column),
#  and indicator of significance level if appropriate (3rd column)
  cat(c("pi0:", qobj$pi0, "\n\n"), file=filename, append=FALSE)
  if(any(names(qobj) == "fdr.level")) {
    cat(c("FDR level:", qobj$fdr.level, "\n\n"), file=filename, append=TRUE)
    cat(c("p-value q-value significant", "\n"), file=filename, append=TRUE) ## change by Alan (space-delimited now)
#    for(i in 1:length(qobj$qval)) {
#      cat(c(qobj$pval[i], "\t", qobj$qval[i], "\t", qobj$significant[i], "\n"), file=filename, append=TRUE)
#    }
    write(t(cbind(qobj$pval, qobj$qval, qobj$significant)), file=filename, ncolumns=3, append=TRUE) ## change by Alan
  }
  else {
    cat(c("p-value q-value", "\n"), file=filename, append=TRUE)
#    for(i in 1:length(qobj$qval)) {
#      cat(c(qobj$pval[i], "\t", qobj$qval[i], "\n"), file=filename, append=TRUE)
#    }
    write(t(cbind(qobj$pval, qobj$qval)), file=filename, ncolumns=2, append=TRUE)
  }
}

qsummary <- function (qobj, cuts=c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1), digits=getOption("digits"), ...) {
  cat("\nCall:\n", deparse(qobj$call), "\n\n", sep = "")
  cat("pi0:",format(qobj$pi0, digits=digits),"\n", sep="\t")
  cat("\n")
  cat("Cumulative number of significant calls:\n")
  cat("\n")
  counts <- sapply(cuts, function(x) c("p-value"=sum(qobj$pvalues < x), "q-value"=sum(qobj$qvalues < x)))
  colnames(counts) <- paste("<", cuts, sep="")
  print(counts)
  cat("\n")
  invisible(qobj)
}

summary.qvalue <- function(object, ...) {
  qsummary(object, ...)
}

qvalue.gui <- function(dummy = NULL) {
    invisible(NULL);
}
