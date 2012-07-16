# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Low-level function for resampling-based prediction error measures
#' 
#' Basic function to estimate the prediction error of a model via (repeated) 
#' \eqn{K}-fold cross-validation, (repeated) random splitting (also known as 
#' random subsampling or Monte Carlo cross-validation), or the bootstrap.  The 
#' model is thereby specified by an unevaluated function call to a model 
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
#' If \code{data} is supplied, all variables required for fitting the models 
#' are added as one argument to the function call, which is the typical 
#' behavior of model fitting functions with a \code{\link[stats]{formula}} 
#' interface.  In this case, a character string specifying the argument name 
#' can be passed via \code{names} (the default is to use \code{"data"}).  
#' 
#' If \code{x} is supplied, on the other hand, the predictor matrix and the 
#' response are added as separate arguments to the function call.  In this 
#' case, \code{names} should be a character vector of length two, with the 
#' first element specifying the argument name for the predictor matrix and the 
#' second element specifying the argument name for the response (the default is 
#' to use \code{c("x", "y")}).  It should be noted that \code{data} takes 
#' precedence over \code{x} if both are supplied.
#' 
#' @param call  an unevaluated function call for fitting a model (see 
#' \code{\link{call}}).
#' @param data  a data frame containing the variables required for fitting the 
#' models.  This is typically used if the model in the function call is 
#' described by a \code{\link[stats]{formula}}.
#' @param x  a numeric matrix containing the predictor variables.  This is 
#' typically used if the function call for fitting the models requires the 
#' predictor matrix and the response to be supplied as separate arguments.
#' @param y  a numeric vector or matrix containing the response.
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
#' 
#' @return If only one replication is requested and the prediction loss 
#' function \code{cost} also returns the standard error, a list is returned, 
#' with the first component containing the estimated prediction errors and the 
#' second component the corresponding estimated standard errors.
#' 
#' Otherwise the return value is a numeric matrix in which each column contains 
#' the respective estimated prediction errors from all replications.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perryTuning}}, 
#' \code{\link{cvFolds}}, \code{\link{foldControl}}, 
#' \code{\link{randomSplits}}, \code{\link{splitControl}}, 
#' \code{\link{bootSamples}}, \code{\link{bootControl}}, 
#' \code{\link{cost}}
#' 
## @example inst/doc/examples/example-perryTool.R
#' 
#' @keywords utilities
#' 
#' @export

perryTool <- function(call, data = NULL, x = NULL, y, cost = rmspe, splits, 
        names = NULL, predictArgs = list(), costArgs = list(),  
        envir = parent.frame()) {
    UseMethod("perryTool", splits)
}

#' @S3method perryTool cvFolds
perryTool.cvFolds <- function(call, data = NULL, x = NULL, y, cost = rmspe, 
        splits, names = NULL, predictArgs = list(), costArgs = list(), 
        envir = parent.frame()) {
    # define an expression that obtains predictions for each test data subset 
    # in one replication of cross-validation
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        cvExpr <- expression(
            lapply(subsets, cvXY, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        cvExpr <- expression(
            lapply(subsets, cvData, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the prediction 
    # and computes the prediction error for one replication of cross-validation
    cvFun <- function(r, keepSE) {
        subsets <- getIndices(splits, r)
        # instead of collecting the results from the folds in the original 
        # order of the observations, the response is re-ordered accordingly
        yHat <- eval(cvExpr)       # obtain predictions from the folds
        yHat <- combineData(yHat)  # combine predictions from the folds
        y <- dataSubset(y, unlist(subsets, use.names=FALSE))  # reorder response
        # compute cost function for predicted values
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform (repeated) cross-validation
    R <- splits$R
    if(R == 1) {
        pe <- cvFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), cvFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)              # combine results
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

#' @S3method perryTool randomSplits
perryTool.randomSplits <- function(call, data = NULL, x = NULL, y, 
        cost = rmspe, splits, names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame()) {
    # define an expression that obtains predictions for the test data in one 
    # replication of random splitting
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        cvExpr <- expression(
            cvXY(subset, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        cvExpr <- expression(
            cvData(subset, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the prediction 
    # and computes the prediction error for one replication of random splitting
    cvFun <- function(r, keepSE) {
        subset <- getIndices(splits, r)
        yHat <- eval(cvExpr)   # obtain predictions for test data
        y <- dataSubset(y, subset)  # response for test data
        # compute cost function for predicted values
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform (repeated) random splitting
    R <- splits$R
    if(R == 1) {
        pe <- cvFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), cvFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)              # combine results
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that take the predictors and the response as separate arguments
cvXY <- function(i, call, x, y, names, predictArgs, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, -i)
    call[[names[2]]] <- dataSubset(y, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(x, i), args=predictArgs)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that have one argument for all the data
cvData <- function(i, call, data, names, predictArgs, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(data, i), args=predictArgs)
}

#' @S3method perryTool bootSamples
perryTool.bootSamples <- function(call, data = NULL, x = NULL, y, 
        cost = rmspe, splits, names = NULL, predictArgs = list(), 
        costArgs = list(), envir = parent.frame()) {
    ## compute out-of-bag prediction error
    # define an expression that obtains the out-of-bag predictions for one 
    # bootstrap sample
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        bootExpr <- expression(
            bootXY(s, call=call, x=x, y=y, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        bootExpr <- expression(
            bootData(s, call=call, data=data, names=names, 
                predictArgs=predictArgs, envir=envir)
        )
    }
    # define a function the evaluates the expression to obtain the out-of-bag 
    # predictions and computes the prediction error for one bootstrap sample 
    bootFun <- function(r, keepSE) {
        s <- getIndices(splits, r)
        yHat <- eval(bootExpr)    # obtain out-of-bag predictions
        y <- dataSubset(y, -s)    # out-of-bag response values
        # compute cost function for obtain out-of-bag predictions
        pe <- computeCost(cost, y, yHat, args=costArgs, keepSE=keepSE)
    }
    # perform bootstrap replications
    R <- splits$R
    if(R == 1 && splits$type != "0.632") {
        pe <- bootFun(R, keepSE=TRUE)
        if(!is.list(pe)) pe <- t(pe)
    } else {
        pe <- lapply(seq_len(R), bootFun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)                # combine results
    }
    ## if requested, compute the 0.632 estimator
    if(splits$type == "0.632") {
        # check if the apparent error has already been computed
        if(is.null(ae <- splits$ae)) {
            # fit the model from all observations
            if(is.null(data)) {
                call[[names[1]]] <- x
                call[[names[2]]] <- y
            } else {
                call[[names]] <- data
            }
            fit <- eval(call, envir)
            # predict the response for all observations
            if(is.null(data)) {
                yHat <- doCall(predict, fit, x, args=predictArgs)
            } else {
                yHat <- doCall(predict, fit, data, args=predictArgs)
            }
            # compute the apparent error with the supplied cost function
            ae <- computeCost(cost, y, yHat, args=costArgs, keepSE=FALSE)
        }
        # compute 0.632 estimator as combination of apparent error and 
        # out-of-bag prediction error
        pe <- sweep(0.632 * pe, 2, 0.368 * ae, "+", check.margin=FALSE)
    }
    # return prediction error (add default names if necessary)
    addNames(pe)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that take the predictors and the 
# response as separate arguments
bootXY <- function(i, call, x, y, names, predictArgs, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, i)
    call[[names[2]]] <- dataSubset(y, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(x, -i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that have one argument for all the data
bootData <- function(i, call, data, names, predictArgs, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(data, -i), args=predictArgs)
}





## compute predictions

# generic function to be extensible
perryPredictions <- function(call, data = NULL, x = NULL, y, 
        names = NULL, splits, predictFun = predict, predictArgs = list(), 
        envir = parent.frame()) {
    UseMethod("perryPredictions", splits)
}

# default method for built-in procedures
perryPredictions.default <- function(call, data = NULL, x = NULL, y, 
        names = NULL, splits, predictFun = predict, predictArgs = list(), 
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
            fun(s, call=call, data=data, x=x, y=y, names=names, 
                predictFun=predictFun, predictArgs=predictArgs, 
                envir=envir)
        })
}

# one replication of cross validation for functions that take the predictors 
# and the response as separate arguments
cvXY <- function(folds, call, data, x, y, names, predictFun, predictArgs, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsXY, call=call, x=x, y=y, names=names, 
        predictFun=predictFun, predictArgs=predictArgs, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# one replication of cross validation for functions that that have one argument 
# for all the data
cvData <- function(folds, call, data, x, y, names, predictFun, predictArgs, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsData, call=call, data=data, names=names, 
        predictFun=predictFun, predictArgs=predictArgs, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that take the predictors and the response as separate arguments
rsXY <- function(i, call, data, x, y, names, predictFun, predictArgs, envir) {
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
rsData <- function(i, call, data, x, y, names, predictFun, predictArgs, envir) {
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
bootXY <- function(i, call, data, x, y, names, predictFun, predictArgs, envir) {
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
bootData <- function(i, call, data, x, y, names, predictFun, predictArgs, envir) {
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
            pe <- computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "randomSplits")) { 
        # function to compute prediction loss for test data
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            pe <- computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "bootSamples")) {
        # function to compute prediction loss for out-of-bag observations
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            pe <- computeCost(cost, dataSubset(y, -s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else stop("invalid data splits")
    # compute prediction loss
    R <- splits$R
    boot0.632 <- inherits(splits, "bootSamples") && splits$type == "0.632"
    if(R == 1 && !boot0.632) {
        # keep standard error if returned by prediction loss function, 
        # otherwise set it to NA
        pe <- fun(r=1, keepSE=TRUE)
        if(!is.list(pe)) pe <- list(pe=pe, se=NA)
    } else {
        # discard standard error if returned by prediction loss function and 
        # recompute it from replications
        pe <- lapply(seq_len(R), fun, keepSE=FALSE)  # for each replication
        pe <- combineData(pe, drop=FALSE)            # combine results
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
            if(keepSE) {
                peNames <- names(pe[[1]])
                pe <- list(sapply(pe, "[[", 1), sapply(pe, "[[", 2))
                names(pe) <- peNames
            } else pe <- sapply(pe, "[[", 1)
        }
    } else {
        pe <- doCall(fun, y, yHat, args=args)
        if(is.list(pe) && !keepSE) pe <- pe[[1]]
    }
    pe
}
