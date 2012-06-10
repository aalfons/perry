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


#' @rdname perryFit
#' @method perryFit function
#' @export

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


#' @rdname perryFit
#' @method perryFit call
#' @export

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
