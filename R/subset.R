# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method subset perry
subset.perry <- function(x, select = NULL, ...) {
    if(is.null(select)) return(x)
    x$pe <- x$pe[select]
    x$se <- x$se[select]
    if(!is.null(reps <- x$reps)) x$reps <- reps[, select, drop=FALSE]
    x
}

#' @S3method subset perrySelect
subset.perrySelect <- function(x, subset = NULL, select = NULL, ...) {
    pe <- x$pe
    se <- x$se
    peNames <- peNames(x)
    reps <- x$reps
    # extract subset of models
    if(is.null(subset)) {
        if(!is.null(select)) {
            if(!is.character(select)) select <- peNames[select]
            x$best <- x$best[select]
            select <- c("Fit", select)  # also select column describing models
            x$pe <- pe[, select, drop=FALSE]
            x$se <- se[, select, drop=FALSE]
            if(!is.null(reps)) x$reps <- reps[, select, drop=FALSE]
        }
    } else {
        if(is.character(subset)) subset <- match(subset, pe$Fit)
        if(inherits(x, "perryTuning")) {
            # extract tuning parameters for the models to keep
            x$tuning <- x$tuning[subset, , drop=FALSE]
        }
        # extract CV results for the models to keep
        if(is.null(select)) {
            pe <- pe[subset, , drop=FALSE]
            se <- se[subset, , drop=FALSE]
        } else {
            if(!is.character(select)) select <- peNames[select]
            select <- c("Fit", select)  # also select column describing models
            pe <- pe[subset, select, drop=FALSE]
            se <- se[subset, select, drop=FALSE]
        }
        fits <- pe$Fit  # models to keep
        haveFactor <- is.factor(fits)
        if(haveFactor) {
            # for a factor, unused levels should be dropped and 
            # remaining levels should be in the right order
            fits <- as.character(fits)
            pe$Fit <- se$Fit <- factor(fits, levels=fits)
        }
        x$pe <- pe
        x$se <- se
        # find best model among the remaining ones
        if(is.null(x$selectBest)) x$selectBest <- "min"
        if(is.null(x$seFactor)) x$seFactor <- NA
        if(nrow(pe) > 0 && ncol(pe) > 1 && !all(is.na(fits))) {
            if(x$selectBest == "min") {
                x$best <- sapply(pe[, -1, drop=FALSE], selectMin)
            } else {
                x$best <- sapply(names(pe)[-1], 
                    function(j) {
                        selectHastie(pe[, j], se[, j], seFactor=x$seFactor)
                    })
            }
        } else x$best <- x$best[integer()]  # this ensures empty integer vector
        # extract the CV replicates for the models to keep
        if(!is.null(reps)) {
            # get list indices of replicates for each model, select the list 
            # components to keep, and flatten the list to an index vector
            indices <- split(seq_len(nrow(reps)), reps$Fit)[subset]
            indices <- unlist(indices, use.names=FALSE)
            # use index vector to extract CV replicates
            if(is.null(select)) {
                reps <- reps[indices, , drop=FALSE]
            } else reps <- reps[indices, select, drop=FALSE]
            if(haveFactor) {
                # for a factor, unused levels should be dropped and 
                # remaining levels should be in the right order
                reps$Fit <- factor(as.character(reps$Fit), levels=fits)
            }
            x$reps <- reps
        }
    }
    x
}
