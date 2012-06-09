# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perryFit <- function(object, ...) UseMethod("perryFit")

#' @S3method perryFit default
perryFit.default <- function(object, data = NULL, x = NULL, y, cost = rmspe, 
        splits = foldControl(), names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame(), seed = NULL, ...) {
    ## extract function call for model fit
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    call <- object$call
    if(is.null(call)) stop("function call for model fitting not available")
    ## for the 0.632 bootstrap estimator, compute the apparent error
    if(inherits(splits, c("bootControl", "bootSamples")) && splits$type == "0.632") {
        n <- nobs(y)
        if(is.null(data)) {
            sx <- "x"
            nx <- nobs(x)
            if(is.null(names)) names <- c("x", "y")
        } else {
            sx <- "data"
            nx <- nobs(data)
            if(is.null(names)) names <- c("data")
        }
        if(!isTRUE(n == nx)) {
            stop(sprintf("'%s' must have %d observations", sx, nx))
        }
        # make sure that .Random.seed exists if no seed is supplied
        if(is.null(seed)) {
            if(!exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) {
                runif(1)
            }
            seed <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
        } else set.seed(seed)
        # compute data splits
        if(is.null(getS3method("perryTool", class(splits), optional=TRUE))) {
            splits <- perrySplits(n, control=splits)
        }
        # predict the response for all observations
        if(is.null(data)) {
            yHat <- doCall(predict, object, x, args=predictArgs)
        } else {
            yHat <- doCall(predict, object, data, args=predictArgs)
        }
        # compute the apparent error with the supplied cost function
        splits$ae <- computeCost(cost, y, yHat, args=costArgs, keepSE=FALSE)
    }
    ## call method for unevaluated function calls
    out <- perryFit(call, data=data, x=x, y=y, cost=cost, splits=splits, 
        names=names, predictArgs=predictArgs, costArgs=costArgs, envir=envir, 
        seed=seed, ...)
    out$call <- matchedCall
    out
}

#' @S3method perryFit function
perryFit.function <- function(object, formula, data = NULL, x = NULL, y, 
        args = list(), cost = rmspe, splits = foldControl(), names = NULL, 
        predictArgs = list(), costArgs = list(), envir = parent.frame(), 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    call <- as.call(c(object, args))  # set up unevaluated function call
    haveFormula <- !missing(formula)
    if(haveFormula || !missing(data)) {
        if(is.null(names)) names <- c("formula", "data")
        if(haveFormula) call[[names[1]]] <- formula
        names <- names[-1]
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        data <- eval(mf, envir)
        if(is.empty.model(attr(data, "terms"))) stop("empty model")
        y <- model.response(data)  # extract response from model frame
    }
    ## call method for unevaluated function calls
    out <- perryFit(call, data=data, x=x, y=y, cost=cost, splits=splits, 
        names=names, predictArgs=predictArgs, costArgs=costArgs, 
        envir=envir, seed=seed)
    out$call <- matchedCall
    out
}

#' @S3method perryFit call
perryFit.call <- function(object, data = NULL, x = NULL, y, cost = rmspe, 
        splits = foldControl(), names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    n <- nobs(y)
    if(is.null(data)) {
        sx <- "x"
        nx <- nobs(x)
    } else {
        sx <- "data"
        nx <- nobs(data)
    }
    if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
    # make sure that .Random.seed exists if no seed is supplied
    if(is.null(seed)) {
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    } else set.seed(seed)
    ## compute data splits
    if(is.null(getS3method("perryTool", class(splits), optional=TRUE))) {
        splits <- perrySplits(n, control=splits)
    }
    ## call workhorse function to estimate the prediction error
    pe <- perryTool(object, data, x, y, cost=cost, splits=splits, names=names, 
        predictArgs=predictArgs, costArgs=costArgs, envir=envir)
    ## compute average results in case of more than one replication and 
    ## prepare standard errors
    if(is.list(pe)) {
        R <- 1  # only one replication
        se <- pe[[2]]
        pe <- pe[[1]]
    } else {
        R <- nrow(pe)
        if(R > 1) {
            reps <- pe
            pe <- apply(reps, 2, mean)
            se <- apply(reps, 2, sd)
        } else {
            pe <- addNames(drop(pe))  # drop() removes column name of 1x1 matrix
            se <- rep.int(NA, length(pe))
            names(se) <- names(pe)
        }
    }
    ## construct return object
    out <- list(splits=splits, pe=pe, se=se)
    if(R > 1) out$reps <- reps
    out$seed <- seed
    out$call <- matchedCall
    class(out) <- "perry"
    out
}
