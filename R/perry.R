# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
perry <- function(object, ...) UseMethod("perry")

## LS regression 
#' @S3method perry lm
perry.lm <- function(object, cost = rmspe, splits = foldControl(), 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    # retrieve data from model fit
    if(is.null(data <- object$model)) {
        haveDataArgument <- !is.null(object$call$data)
        if(haveDataArgument) {
            # try to retrieve data from 'x' and 'y' components
            # this only works if the data argument was used to fit the model
            if(!is.null(x <- object[["x"]]) && !is.null(y <- object$y)) {
                x <- removeIntercept(x)
                data <- data.frame(y, x)
            }
        }
        if(!haveDataArgument || is.null(data)) {
            # try to retrieve data from terms component
            data <- try(model.frame(object$terms), silent=TRUE)
            if(inherits(data, "try-error")) stop("model data not available")
        }
    }
    if(is.null(y <- object$y)) y <- model.response(data)
    ## call function perryFit() to estimate the prediction error
    out <- perryFit(object, data=data, y=y, cost=cost, splits=splits, 
        costArgs=list(...), envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}

## MM and SDMD regression
#' @S3method perry lmrob
perry.lmrob <- function(object, cost = rtmspe, splits = foldControl(), 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    # retrieve data from model fit
    if(is.null(data <- object$model)) {
        haveDataArgument <- !is.null(object$call$data)
        if(haveDataArgument) {
            # try to retrieve data from 'x' and 'y' components
            # this only works if the data argument was used to fit the model
            if(!is.null(x <- object[["x"]]) && !is.null(y <- object$y)) {
                x <- removeIntercept(x)
                data <- data.frame(y, x)
            }
        }
        if(!haveDataArgument || is.null(data)) {
            # try to retrieve data from terms component
            data <- try(model.frame(object$terms), silent=TRUE)
            if(inherits(data, "try-error")) stop("model data not available")
        }
    }
    if(is.null(y <- object$y)) y <- model.response(data)
    ## call function perryFit() to estimate the prediction error
    out <- perryFit(object, data=data, y=y, cost=cost, splits=splits, 
        costArgs=list(...), envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}

## LTS regression
#' @S3method perry lts
perry.lts <- function(object, cost = rtmspe, splits = foldControl(), 
        fit = c("reweighted", "raw", "both"), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    object <- object
    if(is.null(x <- object$X) || is.null(y <- object$Y)) {
        if(is.null(data <- object$model)) {
            if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
            if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
            if(inherits(x, "try-error") || inherits(y, "try-error")) {
                stop("model data not available")
            }
        } else {
            x <- model.matrix(object$terms, data)
            y <- model.response(data)
        }
    }
    # predictor matrix is stored with column for intercept (if any)
    x <- removeIntercept(x)
    ## prepare cross-validation
    # extract function call for model fit
    call <- object$call
    call[[1]] <- as.name("ltsReg")
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    call$intercept <- object$intercept
    ## call function perryFit() to estimate the prediction error
    out <- perryFit(call, x=x, y=y, cost=cost, splits=splits, 
        predictArgs=list(fit=fit), costArgs=list(...), 
        envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


## wrapper functions

## (repeated) K-fold cross validation
#' @export
repCV <- function(object, cost = rmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        grouping = NULL, folds = NULL, ...) {
    ## initializations
    if(is.null(folds)) {
        folds <- foldControl(K, R, type=foldType, grouping=grouping)
    }
    ## call function perry() to estimate the prediction error
    perry(object, cost=cost, splits=folds, ...)
}

## repeated random splitting
#' @export
repRS <- function(object, cost = rmspe, m, R = 1, grouping = NULL, 
        splits = NULL, ...) {
    ## initializations
    if(is.null(splits)) splits <- splitControl(m, R, grouping=grouping)
    ## call function perry() to estimate the prediction error
    perry(object, cost=cost, splits=splits, ...)
}

## bootstrap
#' @export
bootPE <- function(object, cost = rmspe, R = 1, 
        peType = c("0.632", "out-of-bag"), grouping = NULL, 
        samples = NULL, ...) {
    ## initializations
    if(is.null(samples)) {
        samples <- bootControl(R, type=peType, grouping=grouping)
    }
    ## call function perry() to estimate the prediction error
    perry(object, cost=cost, splits=samples, ...)
}
