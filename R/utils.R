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
    } else do.call(rbind, x)
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

# check if a list or object has a certain component
hasComponent <- function(x, name) name %in% names(x)

# check if a generic function has a method for a certain class
# function name needs to be supplied instead of the function itself
hasMethod <- function(fun, class) {
    !is.null(getS3method(fun, class, optional=TRUE))
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
