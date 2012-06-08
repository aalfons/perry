# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method print cvFolds
print.cvFolds <- function(x, ...) {
    # print general information
    cvText <- getPrefix(x)
    if(x$R > 1) {
        cvText <- paste(cvText, "with", x$R, "replications")
    }
    cat(paste("\n", cvText, ":", sep=""))
    # print information on folds (add space between folds and subsets)
    subsets <- x$subsets
    if(x$R == 1) {
        cn <- if(is.null(x$grouping)) "Index" else "Group index"
        nblanks <- 2
    } else {
        cn <- as.character(seq_len(x$R))
        nblanks <- 3
    }
    nblanks <- max(nchar(as.character(subsets[, 1]))-nchar(cn[1]), 0) + nblanks
    cn[1] <- paste(c(rep.int(" ", nblanks), cn[1]), collapse="")
    dimnames(subsets) <- list(Fold=x$which, cn)
    print(subsets, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print randomSplits
print.randomSplits <- function(x, ...) {
    # print general information
    if(x$R == 1) {
        cat("\nRandom split\n")
        cn <- if(is.null(x$grouping)) "Index" else "Group index"
    } else {
        cat(sprintf("\n%d random splits\n", x$R))
        cn <- as.character(seq_len(x$R))
    }
    prefix <- if(is.null(x$grouping)) "Observations" else "Groups"
    cat(sprintf("%s in test data:\n", prefix))
    # print indices of items in test data
    subsets <- x$subsets
    colnames(subsets) <- cn
    print(subsets, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print bootSamples
print.bootSamples <- function(x, ...) {
    # print general information
    if(x$R == 1) {
        postfix <- "sample"
        cn <- if(is.null(x$grouping)) "Index" else "Group index"
    } else {
        postfix <- "samples"
        cn <- as.character(seq_len(x$R))
    }
    prefix <- if(is.null(x$grouping)) "Observations" else "Groups"
    cat(sprintf("\n%s in the bootstrap %s:\n", prefix, postfix))
    # print indices of items in test data
    samples <- x$samples
    colnames(samples) <- cn
    print(samples, ...)
    # return object invisibly
    invisible(x)
}

#' @S3method print perry
#' @S3method print summary.perry
print.perry <- print.summary.perry <- function(x, ...) {
    # print cross-validation results
    cat(getPrefix(x$splits), "results:\n")
    print(x$pe, ...)
    # return object invisibly
    invisible(x)
}


## get prefix for print methods

getPrefix <- function(x) UseMethod("getPrefix")

getPrefix.cvFolds <- function(x) {
    if(x$n == x$K) {
        prefix <- "Leave-one-out CV"
    } else prefix <- sprintf("%d-fold CV", x$K)
}

getPrefix.randomSplits <- function(x) "Monte Carlo CV"

getPrefix.bootSamples <- function(x) {
    if(x$type == "out-of-bag") {
        prefix <- "Out-of-bag bootstrap"
    } else prefix <- sprintf("%s bootstrap", x$type)
}

getPrefix.default <- function(x) "Prediction error"
