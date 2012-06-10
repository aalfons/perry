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

# ----------------------

## utilities for plot functions

# get formula for plot functions
getFormula <- function(left, right, conditional = NULL) {
    if(is.null(conditional)) {
        as.formula(paste(left, "~", right))
    } else as.formula(paste(left, "~", right, "|", conditional))
}


# get data in the correct format for lattice graphics
getLatticeData <- function(x, ...) UseMethod("getLatticeData")

getLatticeData.perry <- function(x, select = NULL, reps = TRUE, 
        seFactor = NA, ...) {
    # extract subset of models
    x <- subset(x, select=select)
    if(reps) {
        PE <- x$reps
        if(is.null(PE)) {
            stop("replications not available")
        } else PE <- as.data.frame(PE)
    } else {
        PE <- as.data.frame(t(x$pe))
    }
    # stack selected results on top of each other
    peName <- defaultNames(1)
    peNames <- peNames(x)
    npe <- npe(x)
    n <- nrow(PE)
    if(npe == 0) {
        # return data frame of NAs if column is not selected
        PE[, peName] <- rep.int(as.numeric(NA), n)
        if(!reps) SE <- as.numeric(NA)
    } else {
        if(!isTRUE(peNames == peName)) {
            PE <- lapply(peNames, 
                function(j) data.frame(Name=rep.int(j, n), PE=PE[, j]))
            PE <- do.call(rbind, PE)
            names(PE) <- c("Name", peName)
        }
        if(!reps) SE <- unname(x$se)
    }
    # return data
    if(reps) {
        PE
    } else {
        halflength <- seFactor * SE
        lower <- PE[, peName] - halflength
        upper <- PE[, peName] + halflength
        list(PE=PE, lower=lower, upper=upper)
    }
}

getLatticeData.perrySelect <- function(x, subset = NULL, select = NULL, 
    reps = TRUE, seFactor = x$seFactor, numericAsFactor = FALSE, ...) {
    # extract subset of models
    x <- subset(x, subset=subset, select=select)
    fits <- fits(x)
    if(reps) {
        PE <- x$reps
        if(is.null(PE)) stop("replications not available")
    } else {
        PE <- x$pe
        SE <- x$se
    }
    # ensure that models are shown in the correct order and drop unused levels
    # ensure that correct values are shown for a numeric tuning parameter
    if(numericAsFactor && is.double(PE$Fit)) {
        PE$Fit <- factor(shingle(PE$Fit), levels=fits)
        if(!reps) SE$Fit <- factor(shingle(SE$Fit), levels=fits)
    } else if(numericAsFactor || !is.numeric(PE$Fit)) {
        PE$Fit <- factor(PE$Fit, levels=fits)
        if(!reps) SE$Fit <- factor(SE$Fit, levels=fits)
    }
    # stack selected results on top of each other
    peName <- defaultNames(1)
    peNames <- peNames(x)
    nfits <- nfits(x)
    npe <- npe(x)
    n <- nrow(PE)
    if(nfits == 0) {
        # no models selected: no column for grouping
        if(isTRUE(peNames == peName) || npe == 0) {
            # return data frame without column for conditional plots and one NA 
            PE <- data.frame(as.numeric(NA))
            names(PE) <- peName
            if(!reps) SE <- as.numeric(NA)
        } else {
            # return data frame with column for conditional plots and NA values
            PE <- data.frame(peNames, rep.int(as.numeric(NA), npe))
            names(PE) <- c("Name", peName)
            if(!reps) SE <- rep.int(as.numeric(NA), npe)
        }
    } else {
        # include column for grouping
        if(npe == 0) {
            # no results selected: no column for conditional plots and NA values
            PE <- PE[, "Fit", drop=FALSE]
            PE[, peName] <- rep.int(as.numeric(NA), n)
            if(!reps) SE <- rep.int(as.numeric(NA), nfits)
        } else {
            # no column for conditional plots if there is only one column of 
            # results with default name
            if(isTRUE(peNames == peName)) {
                if(!reps) SE <- SE[, peName]
            } else {
                PEFit <- PE[, "Fit", drop=FALSE]
                PE <- lapply(peNames, 
                    function(j) cbind(PEFit, Name=rep.int(j, n), PE=PE[, j]))
                PE <- do.call(rbind, PE)
                names(PE) <- c("Fit", "Name", peName)
                if(!reps) SE <- unlist(SE[, peNames], use.names=FALSE)
            }
        }
    }
    # return data
    if(reps) {
        PE
    } else {
        if(is.null(seFactor)) seFactor <- NA
        halflength <- seFactor * SE
        lower <- PE[, peName] - halflength
        upper <- PE[, peName] + halflength
        list(PE=PE, lower=lower, upper=upper)
    }
}

getLatticeData.perryTuning <- function(x, ...) {
    # adjust column specifying the model in case of only one tuning parameter
    if(ncol(x$tuning) == 1) fits(x) <- x$tuning[, 1]
    # call method for class "perrySelect"
    PE <- getLatticeData.perrySelect(x, ...)
    # return data
    PE
}
