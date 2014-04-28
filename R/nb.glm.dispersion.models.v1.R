##' Specify a dispersion model. The parameters of the specified model
##' are to be estimated from the data using the function
##' \code{optim.disp.apl} or \code{optim.disp.pl}.
##'
##' This functions calls \code{disp.fun.<model>} to specify a
##' dispersion model (a list), using output from a call to
##' \code{disp.predictor.<predictor>} as argument list, where
##' \code{<model>} is \code{model} from the input in lower case (one
##' of "nb2", "nbp", "nbq", "nbs" or "step") and \code{<predictor>} is
##' \code{predictor} from the input (one of "pi", "mu", or "rs")
##' 
##' @title (private) Specify a dispersion model
##' @param nb NB data, output from \code{\link{prepare.nb.data}}
##' @param x a matrix, design matrix (specifying the treatment structure).
##' @param model a string giving the name of the disperion model,
##' can be one of "NB2", "NBP", "NBQ", "NBS" or "step" (not case
##' sensitive).
##' @param predictor a string giving the name of the predictor to use
##' in the dispersion model,  can be one of "pi" and "mu", or "rs". 
##' \code{"pi"}, preliminarily estimated mean relative frequencies; \code{"mu"}, preliminarily estimated
##' mean frequencies; \code{"rs"}, row sums.
##' @param subset a list of logical,
##' @param ... additional parameter to \code{disp.fun.*}
##' @return a list, output from the call to the funtion \code{disp.fun.<model>}.
make.disp = function(nb.data, x, model, predictor, subset=filter.mu.pre(nb.data, x), ...) {
    predictor.fun.name = paste("disp.predictor", predictor, sep=".");
    predictor = do.call(predictor.fun.name, list(nb.data=nb.data, x=x));

    disp.fun.name = paste("disp.fun", tolower(model), sep=".");
    do.call(disp.fun.name, c(predictor, list(subset=subset), ...));
}

##' Specify a dispersion model where the parameters of the
##' model will be estimated separately for different groups
##'
##' @title Specify a dispersion model where the parameters of the
##' model will be estimated separately for different groups
##' @param disp.fun 
##' @param grp.ids 
##' @param predictor 
##' @param subset 
##' @param predictor.label 
##' @param ... 
##' @return a list, 
disp.by.group = function(disp.fun, grp.ids, predictor, subset, predictor.label="Predictor", ...) {

    ## grp.ids = rep(1:4, each=3);
    ## grp.ids = rep(1, 3);
    m = length(grp.ids);
    grps = split(1:m, grp.ids);
    n.grps = length(grps);

    res = lapply(grps, function(grp) do.call(disp.fun,
        list(predictor=predictor[,grp,drop=FALSE], subset=subset, predictor.label=predictor.label, ...))); 

    ## This seems tedious
    par.init = res[[1]]$par.init;
    n.pars = length(par.init);
    par.ids  = list(1:n.pars);
    par.lower = res[[1]]$par.lower;
    par.upper = res[[1]]$par.upper;
    funs = list(res[[1]]$fun);
    offset = res[[1]]$offset;

    if (n.grps > 1) {
        for (i in 2:n.grps) {
            par.init = c(par.init, res[[i]]$par.init);
            par.lower = c(par.lower, res[[i]]$par.lower);
            par.upper = c(par.upper, res[[i]]$par.upper);
            funs[[i]] = res[[i]]$fun
            offset = c(offset, res[[i]]$offset);
            par.ids[[i]] = n.pars * (i-1) + 1:length(res[[i]]$par.init);
        }
    }

    env = new.env(parent=baseenv());
    assign("m", nrow(predictor), envir=env)
    assign("n", ncol(predictor), envir=env)
    assign("funs", funs, envir=env)
    assign("grps", grps, envir=env)
    assign("n.grps", n.grps, envir=env)
    assign("par.ids", par.ids, envir=env)
    fun=function(par){
        x = matrix(0, m, n);
        for (i in 1:n.grps) {
            x[,grps[[i]]]=funs[[i]](par[par.ids[[i]]]);
        }
        x
    }
    environment(fun) = env;

    list(name=paste(res[[1]]$name, "by group"),
         fun=fun, par.init=par.init, par.lower=par.lower, par.upper=par.upper, subset=subset,
         predictor=predictor, offset=offset, predictor.label = predictor.label,
         grp.ids = grp.ids)
}

##' @name  Dispersion Models
##'
##' @title (private) Specify a NB2, NBP, NBS, NBS, or STEP dispersion model
##'
##' @details Specify a NBP dispersion model. The parameters of the specified
##' model are to be estimated from the data using the function
##' \code{optim.disp.apl} or \code{optim.disp.pl}.
##'
##' Under the NBP model, the log dispersion is modeled as a linear
##' function of specified predictor with a scalar offset,
##'
##' log(phi) =  par[1] + par[2] * log(predictor/offset).
##'
##' Under this parameterization, par[1] is the dispersion value when
##' the value of predictor equals the offset. This function will
##' return a function (and related settings) to be estimated by either
##' \code{optim.disp.apl} or \code{optim.disp.pl}.  The logical vector
##' \code{subset} specifieds which rows will be used when estimating
##' the paramters (\code{par}) of the dispersion model.
##'
##' Once estimated, the dispersion function will be applied to all
##' values of the \code{predictor} matrix. Care needs to be taken to
##' either avoid \code{NA/Inf} values when preparing the predictor
##' matrix or handle \code{NA/Inf} values afterwards (e.g., when
##' performing hypothesis tests).
##'

##' @rdname disp.fun
##' @param predictor a m-by-n matrix having the same dimensions as the NB counts, predictor of the dispersion. See Details.
##' @param subset a logical vector of length \eqn{m}, specifying the subset of rows to be used when estimating the dispersion model parameters.
##' @param offset a scalar offset.
##' @param predictor.label a string describing the predictor
##' @param par.init a numeric vector, initial values of par.
##' @param label a string character describing the predictor.
##' @param par.lower a numeric vector, lower bounds of the parameter values.
##' @param par.upper a numeric vector, upper bounds of the parameter values.
##' @return a list
##' \item{fun}{a function that takes a vector, \code{par}, as
##' input and outputs a matrix of dispersion values (same dimension as
##' counts)}
##' \item{par.init, par.lower, par.upper}{same as input}
##' \item{subset}{same as input}
##' \item{predictor, offset, predictor.lable}{same as input}
disp.fun.nb2 = function(predictor,
    subset,
    offset=NULL,
    predictor.label="Predictor",
    par.init = -1
    ) {
    ## 2014-04-21
    env = new.env(parent=baseenv());
    assign("d", dim(predictor), envir=env)
    fun=function(par){
        array(exp(par), d);
    }
    environment(fun) = env;

    list(name="NB2",
         fun = fun, par.init = par.init, subset=subset,
         predictor=predictor, predictor.label=predictor.label)
}

##' @rdname disp.fun
disp.fun.nbp = function(predictor,
    subset,
    offset=median(predictor[subset,]),
    predictor.label="Predictor",
    par.init = c(log(0.1), 0),
    par.lower = c(log(1e-20), -1.1),
    par.upper = c(0, 0.1)
    ) {
    ## 2014-04-18
    ## 2014-04-20
    env = new.env(parent=baseenv());
    assign("z", log(predictor/offset), envir=env)
    fun=function(par){
        exp(par[1] + par[2] * z);
    }
    environment(fun) = env;

    list(name="NBP",
         fun = fun, par.init = par.init, par.lower = par.lower, par.upper = par.upper, subset=subset,
         predictor=predictor, offset=offset, predictor.label=predictor.label)
}

##' @rdname disp.fun
disp.fun.nbq = function(predictor,
    subset,
    offset=median(predictor[subset,]),
    predictor.label="Predictor",
    par.init = c(log(0.1), 0, 0),
    par.lower = c(log(1e-20), -1.0, -0.2),
    par.upper = c(0, 1.0, 0.2)
    ) {
    ## 2014-04-18
    ## 2014-04-20
    env = new.env(parent=baseenv());
    assign("z", log(predictor/offset), envir=env)
    fun = function(par) {
        exp(par[1] + par[2] * z + par[3] * z^2);
    }
    environment(fun) = env;

    list(name="NBQ", fun = fun, par.init = par.init, par.lower = par.lower, par.upper = par.upper, subset=subset,
         predictor=predictor, offset=offset, predictor.label =predictor.label)
}

##' @rdname disp.fun
disp.fun.nbs = function(predictor,
    subset,
    offset=NULL,
    predictor.label="Predictor",
    df = 6,
    par.init = rep(-1, df)
    ) {

    ## 2014-04-18
    z = as.vector(log(predictor));

    ## Specify the boundary knots
    Boundary.knots = range(z[subset]);
    ## Boundary.knots = c(1e-4, quantile(z, 0.99));

    ## Specify the knots
    ## It is not clear yet how to select the knots
    zs = sort(z[subset]);
    m = length(zs);
    l = quantile(z[subset], 0.05);
    r = quantile(z[subset], 0.95);

    ##  if (m > 300) {
    ##   l = max(l, zs[101]+0.01);
    ##   r = min(r, zs[m-100]);
    ##  }

    knots = seq(l, r, length=df-2);

    ## Specify the spline basis by the providing the knots
    ## s = ns(z, df = df);
    s = ns(z, knots=knots, Boundary.knots = Boundary.knots, intercept=TRUE);

    ## Specify the dispersion function as a function of par
    ## Create a minimal enviroment for the function
    env = new.env(parent=baseenv());
    assign("s", s, envir=env)
    assign("d", dim(predictor), envir=env)

    fun = function(par) {
        phi = exp(s %*% par);
        dim(phi) = d;
        phi
    }
    environment(fun) = env;

    list(name="NBS", fun = fun, par.init = par.init, subset=subset,
         predictor=predictor, offset=offset, predictor.label =predictor.label)
}

##' @rdname disp.fun
disp.fun.step = function(predictor,
    subset,
    offset=NULL,
    predictor.label="Predictor",
    df = 6,
    knots = NULL,
    par.init = rep(-1, df)
    ) {
    ## 2014-04-21
    z = log(predictor);

    if (is.null(knots)) {
        p = seq(0, 1, length=df+1)[-c(1, df+1)];
        knots = quantile(z[subset], p);
    } else {
        knots = sort(knots);
    }

    ## Specify the starting and ending positions of each step
    l = c(min(z, knots)-0.01, knots);
    r = c(knots, max(z, knots)+0.01);

    ## Identify indices of z values belonging to each step
    ids = list(df);
    m = length(z);
    for (i in 1:df) {
        ids[[i]] = (1:m)[z>=l[i] & z< r[i]];
    }

    d = dim(predictor);

    ## Specify the dispersion function as a function of par
    ## Create a minimal enviroment for the function
    env = new.env(parent=baseenv());
    assign("df", df, envir=env)
    assign("ids", ids, envir=env)
    assign("d", d, envir=env)
    assign("knots", knots, envir=env)

    fun = function(par) {
        lphi = array(NA, d);

        for (i in 1:df) {
            lphi[ids[[i]]] = par[i];
        }

        exp(lphi)

    }
    environment(fun) = env;

    ## Specify the lower and upper bounds for the components of par
    ## Not needed

    list(name="step", fun = fun, par.init = par.init, subset=subset,
         predictor=predictor, predictor.label = predictor.label)
}

##' Create a logical vector specifyfing the subset of rows to be used when estimating the dispersion model
##'
##' @title  Create a logical vector specifyfing the subset of rows to be used when estimating the dispersion model
##' @param nb 
##' @param x 
##' @param mu.lower 
##' @param mu.upper 
##' @return a logical vector specifyfing the subset of rows to be used when estimating the dispersion model
filter.mu.pre = function(nb.data, x, mu.lower=1, mu.upper=Inf, phi.pre=0.1) {
    ## Fit NB regression models to the rows of the count matrix using a
    ## preliminary constant dispersion parameter, phi0.
    obj = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, x, phi=phi.pre);

    ## Obtain preliminary estimates of the mean relative frequencies.
    mu.pre = obj$mu;

    ## Specify the subset of rows to be used when estimating the dispersion model
    subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0;
}


disp.predictor.pi = function(nb.data, x, phi.pre=0.1) {
    ## preliminary constant dispersion parameter, phi0.
    obj = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

    ## Obtain preliminary estimates of the mean relative frequencies.
    mu.pre = obj$mu;
    pi.pre = t(t(mu.pre)/nb.data$eff.lib.sizes);

    ## Replace extremely small pi values 
    ## 11-22-2013
    eps = 1/sum(nb.data$eff.lib.sizes);
    pi.pre[pi.pre<eps] = eps;

    list(predictor = pi.pre, predictor.label = "pi.pre");
}

##' Dispersion precitor
##' @title Dispersion precitor
##' @param nb.data 
##' @param x 
##' @param phi.pre 
##' @param mu.lower 
##' @param mu.upper 
##' @return a logical vector
disp.predictor.mu = function(nb.data, x, phi.pre=0.1, mu.lower=1, mu.upper=Inf) {
    ## Fit NB regression models to the rows of the count matrix using a
    ## preliminary constant dispersion parameter, phi0.
    obj = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

    ## Obtain preliminary estimates of the mean relative frequencies.
    mu.pre = obj$mu;

    ## Replace extremely small pi values 
    ## 11-22-2013
    eps = 1/ncol(nb.data$counts);
    mu.pre[mu.pre<eps] = eps;

    list(predictor = mu.pre, predictor.label = "mu.pre");
}

disp.predictor.rs = function(nb.data, x, phi.pre=0.1, mu.lower=1, mu.upper=Inf) {
    rs = array(rowSums(nb.data$counts), dim(nb.data$counts));
    
    ## Replace extremely small row sums
    eps = 1/ncol(nb.data$counts);
    rs[rs<eps] = eps;

    list(predictor = rs, predictor.label = "row sum");
}
