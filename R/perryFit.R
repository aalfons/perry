# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Resampling-based prediction error for model evaluation
#' 
#' Estimate the prediction error of a model via (repeated) \eqn{K}-fold 
#' cross-validation, (repeated) random splitting (also known as random 
#' subsampling or Monte Carlo cross-validation), or the bootstrap.  It is 
#' thereby possible to supply an object returned by a model fitting function, 
#' a model fitting function itself, or an unevaluated function call to a model 
#' fitting function.
#' 
#' (Repeated) \eqn{K}-fold cross-validation is performed in the following 
#' way.  The data are first split into \eqn{K} previously obtained blocks of 
#' approximately equal size (given by \code{folds}).  Each of the \eqn{K} data 
#' blocks is left out once to fit the model, and predictions are computed for 
#' the observations in the left-out block with the \code{\link[stats]{predict}} 
#' method of the fitted model.  Thus a prediction is obtained for each 
#' observation.  The response variable and the obtained predictions for all 
#' observations are then passed to the prediction loss function \code{cost} to 
#' estimate the prediction error.  For repeated \eqn{K}-fold cross-validation 
#' (as indicated by \code{splits}), this process is replicated and the 
#' estimated prediction errors from all replications are returned.
#' 
#' (Repeated) random splitting is performed similarly.  In each replication, 
#' the data are split into a training set and a test set at random.  Then the 
#' training data is used to fit the model, and predictions are computed for the 
#' test data.  Hence only the response values from the test data and the 
#' corresponding predictions are passed to the prediction loss function 
#' \code{cost}.
#' 
#' For the bootstrap estimator, each bootstrap sample is used as training data 
#' to fit the model.  The out-of-bag estimator uses the observations that do 
#' not enter the bootstrap sample as test data and computes the prediction loss 
#' function \code{cost} for those out-of-bag observations.  The 0.632 estimator 
#' is computed as a linear combination of the out-of-bag estimator and the 
#' prediction loss of the fitted values of the model computed from the full 
#' sample.
#' 
#' In any case, if the response is a vector but the 
#' \code{\link[stats]{predict}} method of the fitted models returns a matrix, 
#' the prediction error is computed for each column.  A typical use case for 
#' this behavior would be if the \code{\link[stats]{predict}} method returns 
#' predictions from an initial model fit and stepwise improvements thereof.
#' 
#' If \code{formula} or \code{data} are supplied, all variables required for 
#' fitting the models are added as one argument to the function call, which is 
#' the typical behavior of model fitting functions with a 
#' \code{\link[stats]{formula}} interface.  In this case, the accepted values 
#' for \code{names} depend on the method.  For the \code{function} method, a 
#' character vector of length two should supplied, with the first element 
#' specifying the argument name for the formula and the second element 
#' specifying the argument name for the data (the default is to use 
#' \code{c("formula", "data")}).  Note that names for both arguments should be 
#' supplied even if only one is actually used.  For the other methods, which do 
#' not have a \code{formula} argument, a character string specifying the 
#' argument name for the data should be supplied (the default is to use 
#' \code{"data"}).  
#' 
#' If \code{x} is supplied, on the other hand, the predictor matrix and the 
#' response are added as separate arguments to the function call.  In this 
#' case, \code{names} should be a character vector of length two, with the 
#' first element specifying the argument name for the predictor matrix and the 
#' second element specifying the argument name for the response (the default is 
#' to use \code{c("x", "y")}).  It should be noted that the \code{formula} or 
#' \code{data} arguments take precedence over \code{x}.
#' 
#' @aliases print.perry
#' 
#' @param object  the fitted model for which to estimate the prediction error, 
#' a function for fitting a model, or an unevaluated function call for fitting 
#' a model (see \code{\link{call}} for the latter).  In the case of a fitted 
#' model, the object is required to contain a component \code{call} that stores 
#' the function call used to fit the model, which is typically the case for 
#' objects returned by model fitting functions.
#' @param formula  a \code{\link[stats]{formula}} describing the model.
#' @param data  a data frame containing the variables required for fitting the 
#' models.  This is typically used if the model in the function call is 
#' described by a \code{\link[stats]{formula}}.
#' @param x  a numeric matrix containing the predictor variables.  This is 
#' typically used if the function call for fitting the models requires the 
#' predictor matrix and the response to be supplied as separate arguments.
#' @param y  a numeric vector or matrix containing the response.
#' @param args  a list of additional arguments to be passed to the model 
#' fitting function.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' (see \code{\link{cost}}).
#' @param splits  an object of class \code{"cvFolds"} (as returned by 
#' \code{\link{cvFolds}}) or a control object of class \code{"foldControl"} 
#' (see \code{\link{foldControl}}) defining the folds of the data for 
#' (repeated) \eqn{K}-fold cross-validation, an object of class 
#' \code{"randomSplits"} (as returned by \code{\link{randomSplits}}) or a 
#' control object of class \code{"splitControl"} (see 
#' \code{\link{splitControl}}) defining random data splits, or an object of 
#' class \code{"bootSamples"} (as returned by \code{\link{bootSamples}}) or a 
#' control object of class \code{"bootControl"} (see \code{\link{bootControl}}) 
#' defining bootstrap samples.
#' @param names  an optional character vector giving names for the arguments 
#' containing the data to be used in the function call (see \dQuote{Details}).
#' @param predictArgs  a list of additional arguments to be passed to the 
#' \code{\link[stats]{predict}} method of the fitted models.
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param envir  the \code{\link{environment}} in which to evaluate the 
#' function call for fitting the models (see \code{\link{eval}}).
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed down.
#' 
#' @returnClass perry
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem pe  a numeric vector containing the respective estimated 
#' prediction errors.  In case of more than one replication, those are average 
#' values over all replications.
#' @returnItem se  a numeric vector containing the respective estimated 
#' standard errors of the prediction loss.
#' @returnItem reps  a numeric matrix in which each column contains the 
#' respective estimated prediction errors from all replications.  This is 
#' only returned in case of more than one replication.
#' @returnItem seed  the seed of the random number generator before estimation 
#' of the prediction error.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryTool}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link{cvFolds}}, 
#' \code{\link{randomSplits}}, \code{\link{bootSamples}}, 
#' \code{\link{cost}}
#' 
## @example inst/doc/examples/example-perryFit.R
#' 
#' @keywords utilities
#' 
#' @export

perryFit <- function(object, ...) UseMethod("perryFit")


#' @rdname perryFit
#' @method perryFit default
#' @export

perryFit.default <- function(object, data = NULL, x = NULL, y, 
        splits = foldControl(), predictFun = predict, predictArgs = list(), 
        cost = rmspe, costArgs = list(), names = NULL, envir = parent.frame(), 
        seed = NULL, ...) {
    ## extract function call for model fit
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    call <- object$call
    if(is.null(call)) stop("function call for model fitting not available")
    ## for the 0.632 bootstrap estimator, compute the apparent error
    if(inherits(splits, c("bootControl", "bootSamples")) && splits$type == "0.632") {
        # predict the response for all observations
        if(is.null(data)) 
            splits$yHat <- doCall(predictFun, object, x, args=predictArgs)
        else splits$yHat <- doCall(predictFun, object, data, args=predictArgs)
    }
    ## call method for unevaluated function calls
    out <- perryFit(call, data=data, x=x, y=y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, cost=cost, 
        costArgs=costArgs, names=names, envir=envir, seed=seed, ...)
    out$call <- matchedCall
    out
}


#' @rdname perryFit
#' @method perryFit function
#' @export

perryFit.function <- function(object, formula, data = NULL, x = NULL, y, 
        args = list(), splits = foldControl(), predictFun = predict, 
        predictArgs = list(), cost = rmspe, costArgs = list(), names = NULL, 
        envir = parent.frame(), seed = NULL, ...) {
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
    out <- perryFit(call, data=data, x=x, y=y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, cost=cost, 
        costArgs=costArgs, names=names, envir=envir, seed=seed, ...)
    out$call <- matchedCall
    out
}


#' @rdname perryFit
#' @method perryFit call
#' @export

perryFit.call <- function(object, data = NULL, x = NULL, y, 
        splits = foldControl(), predictFun = predict, predictArgs = list(), 
        cost = rmspe, costArgs = list(), names = NULL, envir = parent.frame(), 
        seed = NULL, ...) {
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
    if(!is.null(getS3method("perrySplits", class(splits), optional=TRUE))) 
        splits <- perrySplits(n, control=splits)
    ## call workhorse function to compute the predictions
    yHat <- perryPredictions(object, data, x, y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, envir=envir)
    ## call workhorse function to compute the prediction loss
    pe <- perryCost(splits, y, yHat, cost=cost, costArgs=costArgs)
    ## construct return object
    pe <- c(pe, list(splits=splits, y=y, yHat=yHat, 
            seed=seed, call=matchedCall))
    class(pe) <- "perry"
    pe
}


#' @rdname perryFit
#' @method perryFit perry
#' @export

perryFit.perry <- function(object, cost = rmspe, costArgs = list(), ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    ## call workhorse function to compute the prediction loss
    pe <- perryCost(object$splits, object$y, object$yHat, 
        cost=cost, costArgs=costArgs)
    ## construct return object
    object[names(pe)] <- pe
    object$call <- matchedCall
    object
}


## compute predictions

# generic function to be extensible
perryPredictions <- function(call, data = NULL, x = NULL, y, splits, 
    predictFun = predict, predictArgs = list(), names = NULL, 
    envir = parent.frame()) {
    UseMethod("perryPredictions", splits)
}

# default method for built-in procedures
perryPredictions.default <- function(call, data = NULL, x = NULL, y, 
    splits, predictFun = predict, predictArgs = list(), names = NULL, 
    envir = parent.frame()) {
    # define an expression that obtains predictions in one replication
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        if(inherits(splits, "cvFolds")) fun <- cvXY
        else if (inherits(splits, "randomSplits")) fun <- rsXY
        else if(inherits(splits, "bootSamples")) fun <- bootXY
        else stop("invalid data splits")
    } else {
        if(is.null(names)) names <- "data"
        if(inherits(splits, "cvFolds")) fun <- cvData
        else if (inherits(splits, "randomSplits")) fun <- rsData
        else if(inherits(splits, "bootSamples")) fun <- bootData
        else stop("invalid data splits")
    }
    # obtain list of predictions for all replications
    lapply(seq_len(splits$R), function(r) {
            s <- getIndices(splits, r)
            fun(s, call=call, data=data, x=x, y=y, predictFun=predictFun, 
                predictArgs=predictArgs, names=names, envir=envir)
        })
}

# one replication of cross validation for functions that take the predictors 
# and the response as separate arguments
cvXY <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsXY, call=call, x=x, y=y, predictFun=predictFun, 
        predictArgs=predictArgs, names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# one replication of cross validation for functions that that have one argument 
# for all the data
cvData <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsData, call=call, data=data, predictFun=predictFun, 
        predictArgs=predictArgs, names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that take the predictors and the response as separate arguments
rsXY <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, -i)
    call[[names[2]]] <- dataSubset(y, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(x, i), args=predictArgs)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that have one argument for all the data
rsData <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(data, i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that take the predictors and the 
# response as separate arguments
bootXY <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, i)
    call[[names[2]]] <- dataSubset(y, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(x, -i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that have one argument for all the data
bootData <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(data, -i), args=predictArgs)
}


## compute prediction loss

# generic function to be extensible
perryCost <- function(splits, y, yHat, cost, costArgs = list()) {
    UseMethod("perryCost")
}

# default method for built-in procedures
perryCost.default <- function(splits, y, yHat, cost, costArgs = list()) {
    # initializations
    if(inherits(splits, "cvFolds")) {
        # function to compute prediction loss for all observations
        # response needs to be re-oredered according to cross-validation folds
        fun <- function(r, keepSE) {
            s <- unlist(getIndices(splits, r), use.names=FALSE)
            computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "randomSplits")) { 
        # function to compute prediction loss for test data
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "bootSamples")) {
        # function to compute prediction loss for out-of-bag observations
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            computeCost(cost, dataSubset(y, -s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else stop("invalid data splits")
    # compute prediction loss
    R <- splits$R
    boot0.632 <- inherits(splits, "bootSamples") && splits$type == "0.632"
    if(R == 1) {
        # keep standard error if returned by prediction loss function, 
        # otherwise set it to NA
        pe <- fun(r=1, keepSE=!boot0.632)
        if(!is.list(pe)) pe <- list(pe=pe, se=rep.int(NA, length(pe)))
        pe <- addNames(pe)
    } else {
        # discard standard error if returned by prediction loss function and 
        # recompute it from replications
        pe <- lapply(seq_len(R), fun, keepSE=FALSE)  # for each replication
        pe <- addNames(combineData(pe, drop=FALSE))  # combine results
    }
    # if requested, compute the 0.632 bootstrap estimator
    if(boot0.632) {
        # check if the fitted values from the model using all observations 
        # are available
        if(is.null(yHat <- splits$yHat)) 
            stop("fitted values to compute the apparent error are not available")
        else ae <- computeCost(cost, y, yHat, args=costArgs, keepSE=FALSE)
        # compute 0.632 estimator as combination of apparent error and 
        # out-of-bag prediction error
        if(R == 1) pe <- 0.632 * pe + 0.368 * ae
        else pe <- sweep(0.632 * pe, 2, 0.368 * ae, "+", check.margin=FALSE)
    }
    # aggregate results for more than one replication
    if(R > 1) {
        reps <- pe
        pe <- list(pe=apply(reps, 2, mean), se=apply(reps, 2, sd), reps=reps)
    }
    # return list
    pe
}

# compute cost function for predicted values from one replication
computeCost <- function(fun, y, yHat, args = list(), keepSE = TRUE) {
    # if the response is a vector and the predicted values are a matrix, 
    # compute the cost for each column of the matrix of predictions
    if(is.null(dim(y)) && !is.null(dim(yHat))) {
        pe <- apply(yHat, 2, function(x) doCall(fun, y, x, args=args))
        if(is.list(pe)) {
            # cost function returns list of prediction error and standard error
            if(keepSE) 
                pe <- list(pe=sapply(pe, "[[", 1), se=sapply(pe, "[[", 2))
            else pe <- sapply(pe, "[[", 1)
        }
    } else {
        pe <- doCall(fun, y, yHat, args=args)
        if(is.list(pe)) {
            pe <- if(keepSE) list(pe=pe[[1]], se=pe[[2]]) else pe[[1]]
        }
    }
    pe
}
