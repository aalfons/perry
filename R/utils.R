# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## utilities for prediction error functions

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
    if(!check || all(is.na(match(c("Intercept","(Intercept)"), colnames(x))))) {
        cbind("(Intercept)"=rep.int(1, nrow(x)), x)
    } else x
}

# add default names for prediction error results
addNames <- function(x) UseMethod("addNames")
#addNames.list <- function(x) lapply(x, addNames)         # throws error
addNames.list <- function(x) lapply(x, addNames.default)  # workaround
addNames.default <- function(x) {
    if(is.null(p <- ncol(x))) {
        if(is.null(names(x))) names(x) <- defaultNames(length(x))
    } else {
        if(is.null(colnames(x))) colnames(x) <- defaultNames(p)
    }
    x
}

# combine data (used for predictions from PE folds)
combineData <- function(x, drop = TRUE) {
    if(drop && is.null(dim(x[[1]]))) {
        unlist(x)
    } else do.call("rbind", x)
}

# compute cost function for predicted values
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

# retrieve data subsets
dataSubset <- function(x, i, drop = FALSE) {
    if(is.null(dim(x))) {
        x[i]
    } else x[i, , drop=FALSE]
}

# replace data subsets
"dataSubset<-" <- function(x, i, value) {
    if(is.null(dim(x))) {
        x[i] <- value
    } else x[i, ] <- value
    x
}

# default names for prediction error results
defaultNames <- function(p) {
    if(p == 1) {
        "PE"
    } else if(p > 0) {
        paste("PE", seq_len(p), sep="")
    } else character()
}

# default names for model fits
defaultFitNames <- function(m) {
    if(m == 1) {
        "Fit"
    } else if(m > 0) {
        paste("Fit", seq_len(m), sep="")
    } else character()
}

## call a function by either
# 1) simply evaluating a supplied function for the basic arguments if there are
#    no additional arguments in list format
# 2) evaluating a supplied function with 'do.call' if there are additional 
#    arguments in list format
doCall <- function(fun, ..., args = list()) {
    if(length(args) == 0) {
        fun(...)
    } else do.call(fun, c(list(...), args))
}

# retrieve the number of observations
nobs.default <- function(object, ...) {
    n <- nrow(object)                   # matrix or data.frame
    if(is.null(n)) n <- length(object)  # vector
    n
}

## remove intercept column from design matrix
removeIntercept <- function(x, pos) {
    if(missing(pos)) {
        pos <- match(c("Intercept","(Intercept)"), colnames(x), nomatch = 0)
        if(any(pos > 0)) x[, -pos, drop=FALSE] else x
    } else x[, -pos, drop=FALSE]
}

# find which bootstrap samples have all observations in the bag
whichAllInBag <- function(n, samples) {
    indices <- seq_len(n)
    m <- apply(samples, 2, function(i) length(indices[-i]))  # number out-of-bag
    which(m == 0)
}
