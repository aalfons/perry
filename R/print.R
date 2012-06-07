# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method print cvFolds
print.cvFolds <- function(x, ...) {
    # print general information
    if(x$n == x$K) {
        cvText <- "Leave-one-out CV"
    } else {
        cvText <- sprintf("%d-fold CV", x$K)
        if(x$R > 1) {
            cvText <- paste("Repeated", cvText, "with", x$R, "replications")
        }
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
        cat("\nRandom splits\n")
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
